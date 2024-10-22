import numpy as np
from abc import ABC, abstractmethod


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
        return (np.random.random(shape) <= p).astype(np.int32)

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

class Move:
    def __init__(self, array, cost):
        self.array = array
        self.cost = cost
        self.dim = array.shape[0]

    def flip(self):
        self.array *= -1
        self.cost *= -1

    def orientate(self):
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
        return np.all(self.array == 0)

    def copy(self):
        return Move(self.array.copy(), self.cost)

    def __lt__(self, other):
        return self.cost < other.cost

    def __repr__(self):
        return '[' + ','.join(map(str, self.array)) + ']'


def spoly(f, g):
    """Return the s-vector of numpy arrays f and g."""
    s = Move(f.array - g.array, f.cost - g.cost)
    s.orientate()
    return s


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


def is_neg_disjoint(f, g):
    """
    Returns True if the negative support of Moves f and g are disjoint
    """
    return np.all(~((f.array < 0) & (g.array < 0)))


def is_pos_disjoint(f, g):
    """
    Returns True if the positive support of Moves f and g are disjoint
    """
    return np.all(~((f.array > 0) & (g.array > 0)))


def reduce(f, h, negative=False):
    """
    This will use Move f to reduce Move h. This just means producing a new move
    with vector h.array - f.array, and cost h.cost - f.cost.

    If negative is True, it will set the new Move using the sum instead of the difference
    """
    if not negative:
        return Move(h.array - f.array, h.cost - f.cost)
    else:
        return Move(h.array + f.array, h.cost + f.cost)



def reduces(f, h, negative=False):
    """
    This returns True if Move f reduces Move h, which means that f^+ <= h^+.
    If negative is True, then it checks if f^+ <= h^-.
    """
    if not negative:
        return np.all(positive_part(f.array) <= positive_part(h.array))
    else:
        return np.all(positive_part(f.array) <= negative_part(h.array))


def find_reducer(h, F, skip=False, negative=False):
    """
    Given a Move h and a list of Moves F, find a move which reduces h
    or return None if None exists.
    """
    for f in F:
        if skip and np.all(h.array == f.array):
            continue
        if reduces(f, h, negative=negative):
            return f
    return None


def reduce_by_set(g, F, skip=False):
    """Return a remainder of g when reduced by a collection

    Parameters
    ----------
    g : polynomial
        Dividend polynomial.
    F : list
        List of monic divisor polynomials.
    """
    h = g.copy()

    # Reduce Positive
    while True:
        reducer = find_reducer(h, F, skip=skip)
        #print("Reducing by", reducer)
        if reducer is None:
            break
        else:
            h = reduce(h, reducer)
            h.orientate()

    # Reduce Negative
    while True:
        reducer = find_reducer(h, F, skip=skip, negative=True)
        #print("Reducing Negative by", reducer)
        if reducer is None:
            break
        else:
            h = reduce(h, reducer, negative=True)
            h.orientate()

    return h


def update(G, P, f, strategy='none'):
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

    f.orientate()
    G.append(f)
    P.extend(P_)

    return G, P


def minimalize(G):
    """Return a minimal Groebner basis from arbitrary Groebner basis G."""
    Gmin = []
    for f in sorted(G):
        if all(not reduces(f, g) for g in Gmin):
            Gmin.append(f)
    return Gmin


def interreduce(G):
    """Return the reduced Groebner basis from minimal Groebner basis G."""
    Gred = []
    for i in range(len(G)):
        g = reduce_by_set(G[i], G, True)
        Gred.append(g)
    return Gred


def select(G, P, strategy='normal', cost=None):
    """Select and return a pair from P."""
    assert len(G) > 0, "polynomial list must be nonempty"
    assert len(P) > 0, "pair set must be nonempty"

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

    if S is None:
        G = []
        P = []
        for f in F:
            u = Move(f, np.dot(cost, f))
            G, P = update(G, P, u, strategy=elimination)
            print(f'length G = {len(G)}, length P = {len(P)}')
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

    while P:
        i, j = select(G, P, strategy='degree')
        P.remove((i, j))
        s = spoly(G[i], G[j])
        r = reduce_by_set(s, G_)
        #reward = (-1.0 - s['steps']) if rewards == 'additions' else -1.0
        #stats['polynomial_additions'] += s['steps'] + 1
        #stats['total_reward'] += reward
        #stats['discounted_return'] += discount * reward
        #discount *= gamma
        if not r.is_zero():
            G, P = update(G, P, r, strategy=elimination)
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
        print(f'length G = {len(G)}, length P = {len(P)}')

    return interreduce(minimalize(G))

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
        self.cost_dist = cost_dist
        self.x0_dist = x0_dist
        self.elimination = elimination

    def reset(self):
        """Initialize the polynomial list and pair list for a new Groebner basis computation."""
        #if self.sort_input:
            #F.sort(key=lambda f: self.order(f.LM))

        self.G = []                     # the generators in inserted order
        self.G_ = []                    # the reducers if sort_reducers
        self.P = []                     # the pair set

        self.cost = self.cost_dist.sample()
        x0 = self.x0_dist.sample()
        self.x0 = Move(x0, np.dot(x0, self.cost))

        for u in self.basis:
            f = Move(u, np.dot(u, self.cost))
            self.G, self.P = update(self.G, self.P, f, strategy=self.elimination)
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

        self.curr_obj = reduce_by_set(self.x0, self.G_)

        return (self.G, self.P) if self.P else self.reset()

    def step(self, action):
        """Perform one reduction and return the new polynomial list and pair list."""
        i, j = action
        self.P.remove(action)
        s = spoly(self.G[i], self.G[j])
        r = reduce_by_set(s, self.G_)
        if not r.is_zero():
            self.G, self.P = update(self.G, self.P, r, strategy=self.elimination)
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
        reduced = reduce_by_set(self.x0, self.G_)
        return (self.G, self.P), reduced.cost, len(self.P) == 0, {}

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
