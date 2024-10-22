# cython : profile=True
# distutils: language = c++

from libcpp.vector cimport vector
from libcpp.utility cimport pair

import numpy as np
from abc import ABC, abstractmethod


cdef class Move:
    cdef public long[:] array
    cdef public long cost
    cdef public dim

    def __init__(self, long[:] array, long cost):
        self.array = array
        self.cost = cost
        self.dim = array.shape[0]

    def flip(self):
        cdef int i
        for i in range(self.dim):
            self.array[i] *= -1
        self.cost *= -1

    def orientate(self):
        cdef int i

        if self.cost < 0:
            self.flip()
        elif self.cost == 0:
            for i in range(self.dim):
                if self.array[i] > 0:
                    self.flip()
                    return
                elif self.array[i] < 0:
                    return

    def is_zero(self):
        cdef int i
        for i in range(self.dim):
            if self.array[i] != 0:
                return False
        return True

    def copy(self):
        return Move(self.array, self.cost)

    def __lt__(self, other):
        return self.cost < other.cost

    def __repr__(self):
        return '[' + ','.join(map(str, self.array)) + ']'

def orientate(move, cost):
    cdef int i

    if cost < 0:
        return True
    elif cost == 0:
        for i in range(move.shape[0]):
            if move[i] > 0:
                return True
            elif move[i] < 0:
                return False

    return False


cdef class Pairs:
    cdef vector[pair[long,long]] pairs

    def __init__(self):
        pairs = new vector[pair[long,long]]()

    def add_pair(self, long first, long second):
        cdef pair[long,long] a
        a.first = first
        a.second = second
        self.pairs.push_back(a)

    def remove_pair(self, long first, long second):
        cdef pair[long,long] a
        cdef int i
        a.first = first
        a.second = second
        for i in range(self.pairs.size()):
            if (self.pairs[i] == a):
                self.pairs.erase(self.pairs.begin() + i)

    def size(self):
        return self.pairs.size()

    def __getitem__(self, idx):
        return self.pairs[idx]

    def __repr__(self):
        cdef int i;
        s = ''
        for i in range(self.pairs.size()):
            s += f'({self.pairs[i].first},{self.pairs[i].second})\n'
        return s


class Distribution(ABC):
    """
    Abstract base class for a type which generates an object from a distribution
    on demand. Derived classes must implement the following methods:

    sample:
    seed:
    more to be added as needed
    """
    @abstractmethod
    def sample(self):
        pass

    @abstractmethod
    def seed(self, value):
        pass


class NormalDistribution(Distribution):
    def __init__(self, shape, mean, std):
        self.shape = shape
        self.mean = mean
        self.std = std

    def sample(self):
        return np.random.normal(self.mean, self.std, self.shape)

    def seed(self, value):
        np.random.seed(value)


class BernoulliDistribution(Distribution):
    def __init__(self, shape, p):
        self.shape = shape
        self.p = p

    def sample(self):
        return (np.random.random(self.shape) <= self.p).astype(np.int32)

    def seed(self, value):
        np.random.seed(value)


class ConstantDistribution(Distribution):
    """
    Distribution which with singleton support. This will always
    return the same value when sampled.
    """
    def __init__(self, x):
        self.x = x

    def sample(self):
        return self.x

    def seed(self, value):
        pass



def spoly(f, g):
    """Return the s-vector of numpy arrays f and g."""
    return f - g


def LCM(f, g):
    """
    Return the LCM of two monomials represented as np.arrays.
    The arrays represent the exponent vectors of each monomial, and
    the LCM is just given by the maximum
    """
    return np.maximum(f, g)


def positive_part(f):
    """
    Given an array f, this returns the positive part of f,
    which is max(f, 0)
    """
    return np.maximum(f, 0)


def negative_part(f):
    """
    Given an array f, this returns the negative part of f,
    which is max(-f, 0)
    """
    return np.maximum(-f, 0)


cdef is_neg_disjoint(long[:] f, long[:] g):
    """
    Returns True if the negative support of Moves f and g are disjoint
    """
    cdef int i, dim
    dim = f.shape[0]
    for i in range(dim):
        if f[i] < 0 and g[i] < 0:
            return False
    return True

    #return np.all(~((f.array < 0) & (g.array < 0)))


cdef is_pos_disjoint(long[:] f, long[:] g):
    """
    Returns True if the positive support of Moves f and g are disjoint
    """
    cdef int i, dim
    dim = f.shape[0]
    for i in range(dim):
        if f[i] > 0 and g[i] > 0:
            return False
    return True

    #return np.all(~((f.array > 0) & (g.array > 0)))


def reduce(f, h, negative=False):
    """
    This will use Move f to reduce Move h. This just means producing a new move
    with vector h.array - f.array, and cost h.cost - f.cost.

    If negative is True, it will set the new Move using the sum instead of the difference
    """
    if not negative:
        return h - f
    else:
        return h + f



cpdef reduces(long[:] f, long[:] h, bint negative=False):
    """
    This returns True if Move f reduces Move h, which means that f^+ <= h^+.
    If negative is True, then it checks if f^+ <= h^-.
    """
    cdef int length, i
    length = f.shape[0]
    
    if not negative:
        for i in range(length):
            if f[i] > max(h[i], 0):
                return False
        return True
    else:
        for i in range(length):
            if f[i] > max(-h[i], 0):
                return False
        return True


cpdef find_reducer(long[:] h, long[:,:] F, bint skip=False, bint negative=False):
    """
    Given a Move h and a list of Moves F, find a move which reduces h
    or return None if None exists.
    """
    cdef int i, j, dim, n_moves
    n_moves = F.shape[0]
    dim = F.shape[1]

    for i in range(n_moves):
        if skip:
            for j in range(dim):
                if F[i,j] != h[j]:
                    break
            else:
                continue

        if reduces(F[i], h, negative=negative):
            return i
    return None


def reduce_by_set(g, F, g_cost, costs, skip=False, negative=True):
    """Return a remainder of g when reduced by a collection

    Parameters
    ----------
    g : polynomial
        Dividend polynomial.
    F : list
        List of monic divisor polynomials.
    """
    h = g.copy()
    h_cost = g_cost

    # Reduce Positive
    while True:
        reducer_idx = find_reducer(h, F, skip=skip)
        if reducer_idx is None:
            break
        else:
            h -= F[reducer_idx]
            h_cost -= costs[reducer_idx]
            if orientate(h, h_cost):
                h *= -1
                h_cost *= -1

    if negative:
        # Reduce Negative
        while True:
            reducer_idx = find_reducer(h, F, skip=skip, negative=True)
            if reducer_idx is None:
                break
            else:
                h += F[reducer_idx]
                h_cost += costs[reducer_idx]
                if orientate(h, h_cost):
                    h *= -1
                    h_cost *= -1

    return h, h_cost


def update(G, P, costs, f, f_cost, strategy='none'):
    """Return the updated lists of polynomials and pairs when f is added to the basis G.
    The inputs G and P are modified by this function.

    Parameters
    ----------
    G : list
        Current list of polynomial generators.
    P : list
        Current list of s-pairs.
    f : polynomial
        New polynomial to add to the basis.
    strategy : {'gebauermoeller', 'lcm', 'none'}, optional
        Strategy for pair elimination.

        Strategy can be 'none' (eliminate no pairs), 'lcm' (only eliminate pairs that
        fail the LCM criterion), or 'gebauermoeller' (use full Gebauer-Moeller elimination).

    Examples
    --------
    """
    #print("Adding", f)

    if orientate(f, f_cost):
        f *= -1
        f_cost *= -1

    #lmf = f.LM
    #lmG = [g.LM for g in G] if lmG is None else lmG
    m = len(G)

    if strategy == 'none':
        P_ = [(i, m) for i in range(m)]
    elif strategy == 'lcm':
        P_ = [(i, m) for i in range(m) if is_neg_disjoint(G[i], f) and not is_pos_disjoint(G[i], f)]
    elif strategy == 'gebauermoeller':
        raise NotImplementedError
        def can_drop(p):
            i, j = p
            gam = lcm(lmG[i], lmG[j])
            return div(gam, lmf) and gam != lcm(lmG[i], lmf) and gam != lcm(lmG[j], lmf)
        P[:] = [p for p in P if not can_drop(p)]

        lcms = {}
        for i in range(m):
            lcms.setdefault(lcm(lmG[i], lmf), []).append(i)
        min_lcms = []
        P_ = []
        for gam in sorted(lcms.keys(), key=R.order):
            if all(not div(gam, m) for m in min_lcms):
                min_lcms.append(gam)
                if not any(lcm(lmG[i], lmf) == mul(lmG[i], lmf) for i in lcms[gam]):
                    P_.append((lcms[gam][0], m))
        P_.sort(key=lambda p: p[0])
    else:
        raise ValueError('unknown elimination strategy')

    G = np.vstack([G, f])
    costs = np.append(costs, [f_cost])
    # If pairs to add, add them
    if len(P_) > 0:
        P = np.vstack([P, P_])

    return G, P, costs


def minimalize(G, costs):
    """Return a minimal Groebner basis from arbitrary Groebner basis G."""
    Gmin = []
    for f in G:
        if all(not reduces(f, g) for g in Gmin):
            Gmin.append(f)
    return Gmin


def interreduce(G, costs):
    """Return the reduced Groebner basis from minimal Groebner basis G."""
    Gred = []
    for i in range(len(G)):
        g, g_cost = reduce_by_set(G[i], np.array(G), costs[i], costs, True)
        Gred.append(g)
    return Gred


def select(G, P, strategy='normal', cost=None):
    """Select and return a pair from P."""
    assert len(G) > 0, "polynomial list must be nonempty"
    assert P.size() > 0, "pair set must be nonempty"

    if isinstance(strategy, str):
        strategy = [strategy]

    def strategy_key(p, s):
        """Return a sort key for pair p in the strategy s."""
        if s == 'first':
            return p[1], p[0]
        elif s == 'normal':
            lcm = LCM(positive_part(G[p[0]].array), positive_part(G[p[1]].array))
            return np.dot(cost, lcm)
        elif s == 'degree':
            lcm = LCM(positive_part(G[p[0]].array), positive_part(G[p[1]].array))
            return np.sum(lcm)
        elif s == 'random':
            return np.random.rand()
        else:
            raise ValueError('unknown selection strategy')

    return min(P, key=lambda p: tuple(strategy_key(p, s) for s in strategy))


cdef select_degree(long[:,:] G, long[:,:] P):
    cdef int i, j, total, best_i = 0, best_total=1000000
    cdef long dim = G.shape[1]
    cdef int num_pairs = P.shape[0]

    for i in range(num_pairs):
        total = 0
        pi = P[i,0]
        pj = P[i,1]
        for j in range(dim):
            total += max(max(G[pi,j], 0), max(G[pj,j], 0))
        if total < best_total:
            best_total = best_total
            best_i = i
    return best_i


def buchberger(F, cost, S=None, elimination='none'):
    """Return the Groebner basis for the ideal generated by F using Buchberger's algorithm.

    Parameters
    ----------
    F : list
        List of polynomial generators.
    S : list or None, optional
        List of current remaining s-pairs (None indicates no s-pair has been done yet).
    elimination : {'gebauermoeller', 'lcm', 'none'}, optional
        Strategy for pair elimination.
    """

    dim = len(F[0])

    if S is None:
        G = np.empty(shape=(0,dim), dtype=int)
        P = np.empty(shape=(0,2), dtype=int)
        costs = np.empty(shape=(0,), dtype=float)
        for f in F:
            f_cost = np.dot(f, cost)
            G, P, costs = update(G, P, costs, f, f_cost, strategy=elimination)
            #print(f'length G = {len(G)}, length P = {len(P)}')
    else:
        G = F
        P = S

    """
    stats = {'zero_reductions': 0,
             'nonzero_reductions': 0,
             'polynomial_additions': 0,
             'total_reward': 0.0,
             'discounted_return': 0.0}
    discount = 1.0

    if sort_reducers and len(G) > 0:
        order = G[0].ring.order
        G_ = [g.copy() for g in G]
        G_.sort(key=lambda g: order(g.LM))
        lmG_, keysG_ = [g.LM for g in G_], [order(g.LM) for g in G_]
    else:
    """
    G_ = G

    while P.shape[0] > 0:
        action_idx = select_degree(G, P)
        i, j = P[action_idx]
        P = np.vstack([P[:action_idx], P[action_idx+1:]])
        s = spoly(G[i], G[j])
        cost = costs[i] - costs[j]
        if orientate(s, cost):
            s *= -1
            cost *= -1

        r, r_cost = reduce_by_set(s, G_, cost, costs)
        #reward = (-1.0 - s['steps']) if rewards == 'additions' else -1.0
        #stats['polynomial_additions'] += s['steps'] + 1
        #stats['total_reward'] += reward
        #stats['discounted_return'] += discount * reward
        #discount *= gamma
        if not np.all(r == 0):
            G, P, costs = update(G, P, costs, r, r_cost, strategy=elimination)
            #lmG.append(r.LM)
            """
            if sort_reducers:
                key = order(r.LM)
                index = bisect.bisect(keysG_, key)
                G_.insert(index, r.monic())
                lmG_.insert(index, r.LM)
                keysG_.insert(index, key)
            else:
            """
            G_ = G
                #lmG_ = lmG
            #stats['nonzero_reductions'] += 1
        #else:
            #pass
            #stats['zero_reductions'] += 1
        #print(f'length G = {len(G)}, length P = {len(P)}')

    return interreduce(minimalize(G, costs), costs)

class BuchbergerEnv:
    """An environment for computing a Groebner basis using Buchberger's algorithm.

    Parameters
    ----------
    ideal_dist : str or IdealGenerator, optional
        IdealGenerator or string naming the ideal distribution.
    elimination : {'gebauermoeller', 'lcm', 'none'}, optional
        Strategy for pair elimination.

    Examples
    --------
    >>> env = BuchbergerEnv()
    >>> env.seed(123)
    >>> env.reset()
    """

    def __init__(self, initial_basis, cost_dist, x0_dist, elimination='none'):
        self.basis = initial_basis
        self.dim = len(initial_basis[0])
        self.cost_dist = cost_dist
        self.x0_dist = x0_dist
        self.elimination = elimination

    def reset(self):
        """Initialize the polynomial list and pair list for a new Groebner basis computation."""
        #if self.sort_input:
            #F.sort(key=lambda f: self.order(f.LM))

        self.G = np.empty(shape=(0,self.dim), dtype=int)                     # the generators in inserted order
        self.G_ = self.G                    # the reducers if sort_reducers
        self.P = np.empty(shape=(0,2), dtype=int)                # the pair set
        self.costs = np.empty(shape=(0,), dtype=float)

        self.cost_vector = self.cost_dist.sample()
        self.x0 = self.x0_dist.sample()
        self.x0_cost = np.dot(self.x0, self.cost_vector)

        for u in self.basis:
            u_cost = np.dot(u, self.cost_vector)
            self.G, self.P, self.costs = update(self.G, self.P, self.costs, u, u_cost, strategy=self.elimination)
            """
            self.lmG.append(f.LM)
            if self.sort_reducers:
                key = self.order(f.LM)
                index = bisect.bisect(self.keysG_, key)
                self.G_.insert(index, f.monic())
                self.lmG_.insert(index, f.LM)
                self.keysG_.insert(index, key)
            else:
            """
            self.G_ = self.G
                #self.lmG_ = self.lmG

        self.curr_obj = reduce_by_set(self.x0, self.G_, self.x0_cost, self.costs, negative=False)

        return (self.G, self.P)

    def step(self, action_idx):
        """Perform one reduction and return the new polynomial list and pair list."""
        i, j = self.P[action_idx]
        self.P = np.vstack([self.P[:action_idx], self.P[action_idx+1:]])
        s = spoly(self.G[i], self.G[j])
        cost = self.costs[i] - self.costs[j]
        if orientate(s, cost):
            s *= -1
            cost *= -1

        r, r_cost = reduce_by_set(s, self.G_, cost, self.costs)
        if not np.all(r == 0):
            self.G, self.P, self.costs = update(self.G, self.P, self.costs, r, r_cost, strategy=self.elimination)
            """
            self.lmG.append(r.LM)
            if self.sort_reducers:
                key = self.order(r.LM)
                index = bisect.bisect(self.keysG_, key)
                self.G_.insert(index, r.monic())
                self.lmG_.insert(index, r.LM)
                self.keysG_.insert(index, key)
            else:
            """
            self.G_ = self.G
            #self.lmG_ = self.G_
        #reward = -(1.0 + stats['steps']) if self.rewards == 'additions' else -1.0
        reduced, reduced_cost = reduce_by_set(self.x0, self.G_, self.x0_cost, self.costs, negative=False)
        return (self.G, self.P), reduced_cost, self.P.shape[0] == 0, {}

    def seed(self, seed=None):
        self.cost_dist.seed(seed)
        self.x0_dist.seed(seed)

    def value(self, gamma=0.99):
        _, stats = buchberger([g.copy() for g in self.G],
                              S=self.P.copy(),
                              elimination=self.elimination,
                              rewards=self.rewards,
                              sort_reducers=self.sort_reducers,
                              gamma=gamma)
        return stats['discounted_return']
