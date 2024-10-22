std::pair<Monomial, ReduceStats> reduce(const Monomial& g, const std::vector<Binomial>& F) {
  int steps = 0;
  Monomial h = g;

    // Reduce positive
  while (true) {
    bool found_divisor = false;

    for (const Binomial& f : F) {
      if (is_divisible(h, f)) {
                h = reduce_positive(h, f);
                found_divisor = true;
                steps++;
                break;
      }
    }

    if (!found_divisor) {
            // We have reduced positive terms as much as possible
            break;
    }

  }

  return {h, {steps}};
}



std::pair<Binomial, ReduceStats> reduce(const Binomial& g, const std::vector<Binomial>& F, bool skip) {
  int steps = 0;
  Binomial h = g;

    // Reduce positive
  while (true) {
    bool found_divisor = false;

    for (const Binomial& f : F) {
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
            // We have reduced positive terms as much as possible
            break;
    }

  }

    // Reduce negative
  while (true) {
    bool found_divisor = false;

    for (const Binomial& f : F) {
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
            // We have reduced negative terms as much as possible
            break;
    }

  }

  return {h, {steps}};
}


void update(std::vector<Binomial>& G, std::vector<SPair>& P, const Binomial& f, EliminationType elimination) {

  int m = G.size();
  std::vector<SPair> P_;
  
  switch (elimination) {
  case EliminationType::None:
    for (int i = 0; i < m; i++) {
      P_.push_back(SPair{i, m});
    }
    break;
  case EliminationType::LCM:
    for (int i = 0; i < m; i++) {
      if (lcm(G[i].LM(), f.LM()) != G[i].LM() * f.LM())
        P_.push_back(SPair{i, m});
    }
    break;
  case EliminationType::GebauerMoeller:
    auto can_drop = [&G, &f](const SPair& p) {
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



int select_pair(const std::vector<Binomial>& F, const std::vector<SPair>& P, SelectionType selection) {
  std::function<bool(const SPair&, const SPair&)> select;
  bool random = false;
  std::default_random_engine rng;

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
  }

  if (random)
    return rand() % P.size();
  else
  {
    int best_idx = 0;
    SPair best_pair = P[0];
    for (int i=1; i < P.size(); i++)
    {
      if select(P[i], P[idx])
      {
        best_idx = i;
        best_pair = P[i];
      }
    }
    return best_idx;
  }
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


BuchbergerEnv::BuchbergerEnv(std::string markov_filename,
                     std::string x0_filename,
                 EliminationType elimination,
                 bool sort_input,
                 bool sort_reducers,
                     int dimension,
                     std::string costtype,
                     double param_1,
                     double param_2,
                     double param_3)
    : markov_file(markov_filename), x0_file(x0_filename), elimination(elimination),
          sort_input(sort_input), sort_reducers(sort_reducers), dim(dimension),
            cost_type(costtype), param1(param_1), param2(param_2), param3(param_3) {
}

/*
BuchbergerEnv::BuchbergerEnv(const BuchbergerEnv& other)
    : 
      G(other.G), P(other.P), markov_file(other.markov_file), x0_file(other.x0_file), elimination(other.elimination),
      sort_input(other.sort_input), sort_reducers(other.sort_reducers), x0(other.x0),
            curr_cost(other.curr_cost), cost_vector(other.cost_vector), G_(other.G_)  {
}
*/

/*
BuchbergerEnv& BuchbergerEnv::operator=(const BuchbergerEnv& other) {
    // the only non-default copy is ideal_gen
    // ideal_gen = other.ideal_gen->copy();
    G = other.G;
    P = other.P;
    elimination = other.elimination;
    sort_input = other.sort_input;
    sort_reducers = other.sort_reducers;
        x0 = other.x0;
        curr_cost = other.curr_cost;
        cost_vector = other.cost_vector;
    G_ = other.G_;
    return *this;
}
*/

void BuchbergerEnv::read_markov_file() {
    std::ifstream my_file;
    my_file.open(markov_file, std::ifstream::in);
    int rows;
    int nvars;
    my_file >> rows >> nvars;
    dim = nvars;
    for (int i=0; i < rows; i++)
    {
        std::vector<int> row;
        for (int j=0; j < nvars; j++)
        {
            int val;
            my_file >> val;
            row.push_back(val);
        }
        Binomial b(row, cost_vector);
        update(G, P, b, elimination);

    if (sort_reducers)
      G_.insert(std::upper_bound(G_.begin(), G_.end(), b, [](const Binomial& f, const Binomial& g) { return f.LM() < g.LM(); }), b);
    else
      G_.push_back(b);
    }
    my_file.close();
}

void BuchbergerEnv::read_x0_file()
{
    std::ifstream my_file;
    my_file.open(x0_file, std::ifstream::in);
    int nvars;
    my_file >> nvars;

    std::vector<int> vals;
    for (int i = 0; i < nvars; i++)
    {
        int val;
        my_file >> val;
        vals.push_back(val);
    }
    x0 = Monomial(vals, cost_vector);
}

void BuchbergerEnv::generate_cost_vector() {
  if (cost_type == "VertexCover")
    {
        // This gives each vertex variable a normally distributed cost coefficient
        // All other variables have 0 cost associated with them
        // param1 is the number of vertices
        // param2 is the mean of the normal distribution
        // param3 is the standard deviation of the normal distribution

        std::vector<double> cost;
        std::default_random_engine rng;
    std::normal_distribution<double> distribution(param2, param3);

        int i;
        for (i=0; i < param1; i++)
            cost.push_back(distribution(rng));
        for (; i < dim; i++)
            cost.push_back(0);

        cost_vector = cost;
    }
    else if (cost_type == "TotalVariation")
    {
    }
}

void BuchbergerEnv::reset() {
    // TODO: Write what reset means
  G.clear();
  G_.clear();
  P.clear();
    
    // This sets cost_vector
    generate_cost_vector();

    // This sets G, G_, and P
  read_markov_file();
    // Sets x0
    read_x0_file();
    //generate_cost_vector(cost_params);

    // Initialize current observed cost using basis
    auto [x, greedy_stats] = reduce(x0, G_);
    curr_cost = x.get_cost();
}

/*
void BuchbergerEnv::reset() {
  std::vector<Binomial> F = ideal_gen->next();
  if (sort_input)
    std::sort(F.begin(), F.end(), [](const Binomial& f, const Binomial& g) { return f.LM() < g.LM(); });
  G.clear();
  G_.clear();
  P.clear();
  for (const Binomial& f : F) {
    update(G, P, f, elimination);
    if (sort_reducers)
      G_.insert(std::upper_bound(G_.begin(), G_.end(), f, [](const Binomial& f, const Binomial& g) { return f.LM() < g.LM(); }), f);
    else
      G_.push_back(f);
  }
  if (P.empty())
    reset();
}
*/


double BuchbergerEnv::step(SPair action) {
  P.erase(std::remove(P.begin(), P.end(), action), P.end());
  auto [r, stats] = reduce(spoly(G[action.i], G[action.j]), G_, false);
  if (r.is_nonzero()) {
    update(G, P, r, elimination);
    if (sort_reducers)
      G_.insert(std::upper_bound(G_.begin(), G_.end(), r, [](const Binomial& f, const Binomial& g) { return f.LM() < g.LM(); }), r);
    else
      G_.push_back(r);
  }
    // See if there is improvement in the produced output
    auto [x, greedy_stats] = reduce(x0, G_);
    double reward = curr_cost - x.get_cost();
    curr_cost = x.get_cost();

  return reward;
}


double BuchbergerEnv::value(std::string strategy, double gamma) const {
  if (strategy == "sample") {
    auto [G_, stats] = buchberger(G, P, SelectionType::Degree, elimination, sort_reducers, gamma);
    double best = stats.discounted_return;
    for (int i = 0; i < 100; i++) {
      std::tie(G_, stats) = buchberger(G, P, SelectionType::Random, elimination, sort_reducers, gamma);
      best = std::max(best, stats.discounted_return);
    }
    return best;
  }
  std::map<std::string, SelectionType> select = {
      {"first", SelectionType::First},
      {"degree", SelectionType::Degree},
      {"normal", SelectionType::Normal},
      {"sugar", SelectionType::Sugar},
      {"random", SelectionType::Random}
  };
  auto [G_, stats] = buchberger(G, P, select[strategy], elimination, sort_reducers, gamma);
  return stats.discounted_return;
}


std::vector<int> lead_monomials_vector(const Binomial& f) {
  std::vector<int> lead;
    for (int j = 0; j < f.size(); j++) {
        lead.push_back(f[j]);
    }
  return lead;
}


LeadMonomialsEnv::LeadMonomialsEnv(std::string markov_file,
                     std::string x0_file,
                     bool sort_input,
                     bool sort_reducers,
                     int dimension,
                     std::string costtype,
                     double param_1,
                     double param_2,
                     double param_3)
  : env{markov_file, x0_file, EliminationType::GebauerMoeller, sort_input, sort_reducers, dimension, costtype, param_1, param_2, param_3}
{
  n = env.nvars();
  cols = 2 * env.nvars();
}


void LeadMonomialsEnv::reset() {
  env.reset();
  state.clear();
  leads.clear();
  basis.clear();
  pairs.clear();
  for (const Binomial& g : env.G) {
    std::vector<int> lmv = lead_monomials_vector(g);
    leads.push_back(lmv);
    basis.insert(basis.end(), lmv.begin(), lmv.end());
  }
  for (const auto& p : env.P) {
    state.insert(state.end(), leads[p.i].begin(), leads[p.i].end());
    state.insert(state.end(), leads[p.j].begin(), leads[p.j].end());
        pairs.push_back(p.i);
        pairs.push_back(p.j);
  }
}


double LeadMonomialsEnv::step(int action) {
  double reward = env.step(env.P[action]);
  if (leads.size() < env.G.size())
  {
    std::vector<int> lmv = lead_monomials_vector(env.G[env.G.size()-1]);
    leads.push_back(lmv);
    basis.insert(basis.end(), lmv.begin(), lmv.end());
  }
  state.clear();
  pairs.clear();
  for (const auto& p : env.P) {
    state.insert(state.end(), leads[p.i].begin(), leads[p.i].end());
    state.insert(state.end(), leads[p.j].begin(), leads[p.j].end());
    pairs.push_back(p.i);
    pairs.push_back(p.j);
  }  
  return reward;
}
