
def greedy_algorithm(point, moves):
    steps = 0;
    x = point.copy();

    while True:
        for move in moves:
            if is_divisible(x, f):
                x = reduce_positive(x, f)
                found_divisor = true
                steps++
        else:
            break
  return h, steps

def reduce(g, moves, skip) {
    steps = 0
    u = g.copy()

    while True:
        for move in moves:
            if (skip && np.all(g == move))
                continue

            if (is_divisible(u, f)) {
                u = reduce_positive(u, f)
                found_divisor = true
                steps++
                break
        else:
            break

    while True:
        for move in moves:
            if (skip && np.all(g == move))
                continue

            if (is_divisible_negative(u, f)) {
                u = reduce_negative(u, f)
                steps++
                break

        else:
            break
  return u, steps


def update(moves, pairs, move, elimination) {

    m = len(moves)
    P_ = []
  
    if elimination == 'None':
        for i in range(m): 
            P_.append((i, m))
    elif elimination == 'LCM'
        for i in range(m): {
            if (lcm(pos_part(G[i]), pos_part(f)) != pos_part(G[i]) * pos_part(f)):
                P_.append((i, m))
    elif elimination == 'GebauerMoller':
        def can_drop(i, j): = [&G, &f](const SPair& p) {
              Monomial m = lcm(G[p.i].LM(), G[p.j].LM());
              return (is_divisible(m, f.LM()) &&
                  (m != lcm(G[p.i].LM(), f.LM())) &&
                  (m != lcm(G[p.j].LM(), f.LM())));
                };
    P.erase(std::remove_if(P.begin(), P.end(), can_drop), P.end());

    std::map<Monomial, std::vector<int>> lcms;
    for (int i = 0; i < m; i++) {
      lcms[lcm(G[i].LM(), f.LM())].push_back(i);
    }
    std::vector<Monomial> min_lcms;
    for (const auto& pair : lcms) {  // will be in sorted order because std::map
      Monomial mon = pair.first;
      std::vector<int> v = pair.second;
      if (std::all_of(min_lcms.begin(), min_lcms.end(), [&mon](const Monomial& m) { return !is_divisible(mon, m); })) {
    min_lcms.push_back(mon);
    if (std::none_of(v.begin(), v.end(), [&G, &f](int i) { return lcm(G[i].LM(), f.LM()) == G[i].LM() * f.LM(); }))
      P_.push_back(SPair{v[0], m});
      }
    }
    std::sort(P_.begin(), P_.end(), [](const SPair& p1, const SPair& p2) { return p1.i < p2.i; });

    break;
  }

  G.push_back(f);
  P.insert(P.end(), P_.begin(), P_.end());
}


std::vector<Binomial> minimalize(const std::vector<Binomial>& G) {
  std::vector<Binomial> G_ = G;
  std::sort(G_.begin(), G_.end(), [](const Binomial& f, const Binomial& g) { return f.LM() < g.LM(); });
  std::vector<Binomial> Gmin;
  for (const Binomial& g : G_) {
    if (std::none_of(Gmin.begin(), Gmin.end(), [&g](const Binomial& f) { return is_divisible(g.LM(), f.LM()); })) {
      Gmin.push_back(g);
    }
  }
  return Gmin;
}


std::vector<Binomial> interreduce(const std::vector<Binomial>& G) {
  std::vector<Binomial> Gred;
  for (const Binomial& g : G) {
        // We reduce g by G, we need to skip over the element g
    Gred.push_back(reduce(g, G, true).first);
  }
  return Gred;
}


std::pair<std::vector<Binomial>, BuchbergerStats> buchberger(const std::vector<Binomial>& F,
                                   SelectionType selection,
                                   EliminationType elimination,
                                   bool sort_input,
                                   bool sort_reducers,
                                   double gamma,
                                   std::optional<int> seed) {
  std::vector<Binomial> G;
  std::vector<SPair> P;
  for (const Binomial& f : F) {
    update(G, P, f, elimination);
  }

  return buchberger(G, P, selection, elimination, sort_reducers, gamma, seed);
}


std::pair<std::vector<Binomial>, BuchbergerStats> buchberger(const std::vector<Binomial>& F,
                                   const std::vector<SPair>& S,
                                   SelectionType selection,
                                   EliminationType elimination,
                                   bool sort_reducers,
                                   double gamma,
                                   std::optional<int> seed) {
  std::vector<Binomial> G = F;
  std::vector<Binomial> G_ = F;
  std::vector<SPair> P = S;
  BuchbergerStats stats = {};
  double discount = 1.0;

  if (sort_reducers)
    std::sort(G_.begin(), G_.end(), [](const Binomial& f, const Binomial& g) { return f.LM() < g.LM(); });

  std::function<bool(const SPair&, const SPair&)> select;
  bool random = false;
  std::default_random_engine rng;
  std::random_device rand;
  switch (selection) {
  case SelectionType::First:
    select = [](const SPair& p1, const SPair& p2) {
           return std::tie(p1.j, p1.i) < std::tie(p2.j, p2.i);
         };
    random = false;
    break;
  case SelectionType::Degree:
    select = [&G](const SPair& p1, const SPair& p2) {
           int d1 = lcm(G[p1.i].LM(), G[p1.j].LM()).deg();
           int d2 = lcm(G[p2.i].LM(), G[p2.j].LM()).deg();
           return std::tie(d1, p1.j, p1.i) < std::tie(d2, p2.j, p2.i);
         };
    random = false;
    break;
  case SelectionType::Normal:
    select = [&G](const SPair& p1, const SPair& p2) {
           Monomial m1 = lcm(G[p1.i].LM(), G[p1.j].LM());
           Monomial m2 = lcm(G[p2.i].LM(), G[p2.j].LM());
           return std::tie(m1, p1.j, p1.i) < std::tie(m2, p2.j, p2.i);
         };
    random = false;
    break;
  case SelectionType::Sugar:
    select = [&G](const SPair& p1, const SPair& p2) {
           Monomial m1 = lcm(G[p1.i].LM(), G[p1.j].LM());
           Monomial m2 = lcm(G[p2.i].LM(), G[p2.j].LM());
           int s1 = std::max(G[p1.i].sugar() + (m1 / G[p1.i].LM()).deg(),
                 G[p1.j].sugar() + (m1 / G[p1.j].LM()).deg());
           int s2 = std::max(G[p2.i].sugar() + (m2 / G[p2.i].LM()).deg(),
                 G[p2.j].sugar() + (m2 / G[p2.j].LM()).deg());
           return std::tie(s1, m1, p1.j, p1.i) < std::tie(s2, m2, p2.j, p2.i);
         };
    random = false;
    break;
  case SelectionType::Random:
    if (seed) {
      rng.seed(seed.value());
    } else {
      rng.seed(rand());
    }
    random = true;
    break;
  case SelectionType::Last:
    select = [](const SPair& p1, const SPair& p2) {
           return std::tie(p1.j, p1.i) > std::tie(p2.j, p2.i);
         };
    random = false;
    break;
  case SelectionType::Codegree:
    select = [&G](const SPair& p1, const SPair& p2) {
           int d1 = lcm(G[p1.i].LM(), G[p1.j].LM()).deg();
           int d2 = lcm(G[p2.i].LM(), G[p2.j].LM()).deg();
           return std::tie(d1, p1.j, p1.i) > std::tie(d2, p2.j, p2.i);
         };
    random = false;
    break;
  case SelectionType::Strange:
    select = [&G](const SPair& p1, const SPair& p2) {
           Monomial m1 = lcm(G[p1.i].LM(), G[p1.j].LM());
           Monomial m2 = lcm(G[p2.i].LM(), G[p2.j].LM());
           return std::tie(m1, p1.j, p1.i) > std::tie(m2, p2.j, p2.i);
         };
    random = false;
    break;
  case SelectionType::Spice:
    select = [&G](const SPair& p1, const SPair& p2) {
           Monomial m1 = lcm(G[p1.i].LM(), G[p1.j].LM());
           Monomial m2 = lcm(G[p2.i].LM(), G[p2.j].LM());
           int s1 = std::max(G[p1.i].sugar() + (m1 / G[p1.i].LM()).deg(),
                 G[p1.j].sugar() + (m1 / G[p1.j].LM()).deg());
           int s2 = std::max(G[p2.i].sugar() + (m2 / G[p2.i].LM()).deg(),
                 G[p2.j].sugar() + (m2 / G[p2.j].LM()).deg());
           return std::tie(s1, m1, p1.j, p1.i) > std::tie(s2, m2, p2.j, p2.i);
         };
    random = false;
    break;
  }

  while (!P.empty()) {
    auto iter = random ? choice(P.begin(), P.end(), rng) : std::min_element(P.begin(), P.end(), select);
    SPair p = *iter;
    P.erase(iter);
    auto [r, s] = reduce(spoly(G[p.i], G[p.j]), G_, false);
        // TODO: Update Reward criterion here
    double reward = (-1.0 - s.steps);
    stats.polynomial_additions += s.steps + 1;
    stats.total_reward += reward;
    stats.discounted_return += discount * reward;
    discount *= gamma;
    if (r.is_nonzero()) {
      update(G, P, r, elimination);
      stats.nonzero_reductions++;
      if (sort_reducers)
                G_.insert(std::upper_bound(G_.begin(), G_.end(), r, [](const Binomial& f, const Binomial& g) { return f.LM() < g.LM(); }), r);
      else
                G_.push_back(r);
    } else {
      stats.zero_reductions++;
    }
  }

  return {interreduce(minimalize(G)), stats};
}
/*
 * an environment for computing groebner bases with buchberger's algorithm.
 */

#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <vector>
#include <fstream>
#include <iostream>

#include "buchberger.h"
//#include "ideals.h"
#include "polynomials.h"


std::pair<monomial, reducestats> reduce(const monomial& g, const std::vector<binomial>& f) {
  int steps = 0;
  monomial h = g;

    // reduce positive
  while (true) {
    bool found_divisor = false;

    for (const binomial& f : f) {
      if (is_divisible(h, f)) {
                h = reduce_positive(h, f);
                found_divisor = true;
                steps++;
                break;
      }
    }

    if (!found_divisor) {
            // we have reduced positive terms as much as possible
            break;
    }

  }

  return {h, {steps}};
}



std::pair<binomial, reducestats> reduce(const binomial& g, const std::vector<binomial>& f, bool skip) {
  int steps = 0;
  binomial h = g;

    // reduce positive
  while (true) {
    bool found_divisor = false;

    for (const binomial& f : f) {
            if (skip && (g == f))
                continue;

      if (is_divisible(h, f)) {
                h = reduce_positive(h, f);
                found_divisor = true;
                steps++;
                break;
      }
    }

    if (!found_divisor) {
            // we have reduced positive terms as much as possible
            break;
    }

  }

    // reduce negative
  while (true) {
    bool found_divisor = false;

    for (const binomial& f : f) {
            if (skip && (g == f))
                continue;

      if (is_divisible_negative(h, f)) {
                h = reduce_negative(h, f);
                found_divisor = true;
                steps++;
                break;
      }
    }

    if (!found_divisor) {
            // we have reduced negative terms as much as possible
            break;
    }

  }

  return {h, {steps}};
}


void update(std::vector<binomial>& g, std::vector<spair>& p, const binomial& f, eliminationtype elimination) {

  int m = g.size();
  std::vector<spair> p_;
  
  switch (elimination) {
  case eliminationtype::none:
    for (int i = 0; i < m; i++) {
      p_.push_back(spair{i, m});
    }
    break;
  case eliminationtype::lcm:
    for (int i = 0; i < m; i++) {
      if (lcm(g[i].lm(), f.lm()) != g[i].lm() * f.lm())
                p_.push_back(spair{i, m});
    }
    break;
  case eliminationtype::gebauermoeller:
    auto can_drop = [&g, &f](const spair& p) {
              monomial m = lcm(g[p.i].lm(), g[p.j].lm());
              return (is_divisible(m, f.lm()) &&
                  (m != lcm(g[p.i].lm(), f.lm())) &&
                  (m != lcm(g[p.j].lm(), f.lm())));
                };
    p.erase(std::remove_if(p.begin(), p.end(), can_drop), p.end());

    std::map<monomial, std::vector<int>> lcms;
    for (int i = 0; i < m; i++) {
      lcms[lcm(g[i].lm(), f.lm())].push_back(i);
    }
    std::vector<monomial> min_lcms;
    for (const auto& pair : lcms) {  // will be in sorted order because std::map
      monomial mon = pair.first;
      std::vector<int> v = pair.second;
      if (std::all_of(min_lcms.begin(), min_lcms.end(), [&mon](const monomial& m) { return !is_divisible(mon, m); })) {
    min_lcms.push_back(mon);
    if (std::none_of(v.begin(), v.end(), [&g, &f](int i) { return lcm(g[i].lm(), f.lm()) == g[i].lm() * f.lm(); }))
      p_.push_back(spair{v[0], m});
      }
    }
    std::sort(p_.begin(), p_.end(), [](const spair& p1, const spair& p2) { return p1.i < p2.i; });

    break;
  }

  g.push_back(f);
  p.insert(p.end(), p_.begin(), p_.end());
}


std::vector<binomial> minimalize(const std::vector<binomial>& g) {
  std::vector<binomial> g_ = g;
  std::sort(g_.begin(), g_.end(), [](const binomial& f, const binomial& g) { return f.lm() < g.lm(); });
  std::vector<binomial> gmin;
  for (const binomial& g : g_) {
    if (std::none_of(gmin.begin(), gmin.end(), [&g](const binomial& f) { return is_divisible(g.lm(), f.lm()); })) {
      gmin.push_back(g);
    }
  }
  return gmin;
}


std::vector<binomial> interreduce(const std::vector<binomial>& g) {
  std::vector<binomial> gred;
  for (const binomial& g : g) {
        // we reduce g by g, we need to skip over the element g
    gred.push_back(reduce(g, g, true).first);
  }
  return gred;
}


std::pair<std::vector<binomial>, buchbergerstats> buchberger(const std::vector<binomial>& f,
                                   selectiontype selection,
                                   eliminationtype elimination,
                                   bool sort_input,
                                   bool sort_reducers,
                                   double gamma,
                                   std::optional<int> seed) {
  std::vector<binomial> g;
  std::vector<spair> p;
  for (const binomial& f : f) {
    update(g, p, f, elimination);
  }

  return buchberger(g, p, selection, elimination, sort_reducers, gamma, seed);
}


std::pair<std::vector<binomial>, buchbergerstats> buchberger(const std::vector<binomial>& f,
                                   const std::vector<spair>& s,
                                   selectiontype selection,
                                   eliminationtype elimination,
                                   bool sort_reducers,
                                   double gamma,
                                   std::optional<int> seed) {
  std::vector<binomial> g = f;
  std::vector<binomial> g_ = f;
  std::vector<spair> p = s;
  buchbergerstats stats = {};
  double discount = 1.0;

  if (sort_reducers)
    std::sort(g_.begin(), g_.end(), [](const binomial& f, const binomial& g) { return f.lm() < g.lm(); });

  std::function<bool(const spair&, const spair&)> select;
  bool random = false;
  std::default_random_engine rng;
  std::random_device rand;
  switch (selection) {
  case selectiontype::first:
    select = [](const spair& p1, const spair& p2) {
           return std::tie(p1.j, p1.i) < std::tie(p2.j, p2.i);
         };
    random = false;
    break;
  case selectiontype::degree:
    select = [&g](const spair& p1, const spair& p2) {
           int d1 = lcm(g[p1.i].lm(), g[p1.j].lm()).deg();
           int d2 = lcm(g[p2.i].lm(), g[p2.j].lm()).deg();
           return std::tie(d1, p1.j, p1.i) < std::tie(d2, p2.j, p2.i);
         };
    random = false;
    break;
  case selectiontype::normal:
    select = [&g](const spair& p1, const spair& p2) {
           monomial m1 = lcm(g[p1.i].lm(), g[p1.j].lm());
           monomial m2 = lcm(g[p2.i].lm(), g[p2.j].lm());
           return std::tie(m1, p1.j, p1.i) < std::tie(m2, p2.j, p2.i);
         };
    random = false;
    break;
  case selectiontype::sugar:
    select = [&g](const spair& p1, const spair& p2) {
           monomial m1 = lcm(g[p1.i].lm(), g[p1.j].lm());
           monomial m2 = lcm(g[p2.i].lm(), g[p2.j].lm());
           int s1 = std::max(g[p1.i].sugar() + (m1 / g[p1.i].lm()).deg(),
                 g[p1.j].sugar() + (m1 / g[p1.j].lm()).deg());
           int s2 = std::max(g[p2.i].sugar() + (m2 / g[p2.i].lm()).deg(),
                 g[p2.j].sugar() + (m2 / g[p2.j].lm()).deg());
           return std::tie(s1, m1, p1.j, p1.i) < std::tie(s2, m2, p2.j, p2.i);
         };
    random = false;
    break;
  case selectiontype::random:
    if (seed) {
      rng.seed(seed.value());
    } else {
      rng.seed(rand());
    }
    random = true;
    break;
  case selectiontype::last:
    select = [](const spair& p1, const spair& p2) {
           return std::tie(p1.j, p1.i) > std::tie(p2.j, p2.i);
         };
    random = false;
    break;
  case selectiontype::codegree:
    select = [&g](const spair& p1, const spair& p2) {
           int d1 = lcm(g[p1.i].lm(), g[p1.j].lm()).deg();
           int d2 = lcm(g[p2.i].lm(), g[p2.j].lm()).deg();
           return std::tie(d1, p1.j, p1.i) > std::tie(d2, p2.j, p2.i);
         };
    random = false;
    break;
  case selectiontype::strange:
    select = [&g](const spair& p1, const spair& p2) {
           monomial m1 = lcm(g[p1.i].lm(), g[p1.j].lm());
           monomial m2 = lcm(g[p2.i].lm(), g[p2.j].lm());
           return std::tie(m1, p1.j, p1.i) > std::tie(m2, p2.j, p2.i);
         };
    random = false;
    break;
  case selectiontype::spice:
    select = [&g](const spair& p1, const spair& p2) {
           monomial m1 = lcm(g[p1.i].lm(), g[p1.j].lm());
           monomial m2 = lcm(g[p2.i].lm(), g[p2.j].lm());
           int s1 = std::max(g[p1.i].sugar() + (m1 / g[p1.i].lm()).deg(),
                 g[p1.j].sugar() + (m1 / g[p1.j].lm()).deg());
           int s2 = std::max(g[p2.i].sugar() + (m2 / g[p2.i].lm()).deg(),
                 g[p2.j].sugar() + (m2 / g[p2.j].lm()).deg());
           return std::tie(s1, m1, p1.j, p1.i) > std::tie(s2, m2, p2.j, p2.i);
         };
    random = false;
    break;
  }

  while (!p.empty()) {
    auto iter = random ? choice(p.begin(), p.end(), rng) : std::min_element(p.begin(), p.end(), select);
    spair p = *iter;
    p.erase(iter);
    auto [r, s] = reduce(spoly(g[p.i], g[p.j]), g_, false);
        // todo: update reward criterion here
    double reward = (-1.0 - s.steps);
    stats.polynomial_additions += s.steps + 1;
    stats.total_reward += reward;
    stats.discounted_return += discount * reward;
    discount *= gamma;
    if (r.is_nonzero()) {
      update(g, p, r, elimination);
      stats.nonzero_reductions++;
      if (sort_reducers)
                g_.insert(std::upper_bound(g_.begin(), g_.end(), r, [](const binomial& f, const binomial& g) { return f.lm() < g.lm(); }), r);
      else
                g_.push_back(r);
    } else {
      stats.zero_reductions++;
    }
  }

  return {interreduce(minimalize(g)), stats};
}
