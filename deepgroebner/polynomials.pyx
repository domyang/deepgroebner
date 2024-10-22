import numpy as np
cimport numpy as np
np.import_array()

INT_TYPE = np.int64
DOUBLE_TYPE = np.float64
ctypedef np.int64_t INT_TYPE_t
ctypedef np.float64_t DOUBLE_TYPE_t


cdef class Binomial:
    cdef public INT_TYPE_t[:] exponent
    cdef public int degree
    cdef public double cost
    cdef public DOUBLE_TYPE_t[:] cost_vector

cdef class Monomial:
    
    cdef public INT_TYPE_t[:] exponent
    cdef public int degree
    cdef public double cost
    cdef public DOUBLE_TYPE_t[:] cost_vector

    def __init__(self, INT_TYPE_t[:] exponent, DOUBLE_TYPE_t[:] cost_vector):
        self.exponent = exponent
        self.degree = exponent.shape[0]
        self.cost = np.dot(exponent, cost_vector)
        self.cost_vector = cost_vector

    def __getitem__(self, int idx):
        return self.exponent[idx]

    @property
    def size(self):
        return self.exponent.shape[0]

    def __mul__(self, other):
        exp_sum = np.zeros(self.size, dtype=INT_TYPE)
        cdef INT_TYPE_t[:] exp_sum_view = exp_sum
        cdef int i
        for i in range(self.size):
            exp_sum_view[i] = self.exponent[i] + other.exponent[i]
        return Monomial(exp_sum, self.cost_vector)

    def __div__(self, other):
        exp_diff = np.zeros(self.size, dtype=INT_TYPE)
        cdef INT_TYPE_t[:] exp_diff_view = exp_diff
        cdef int i
        for i in range(self.size):
            exp_diff_view[i] = self.exponent[i] - other.exponent[i]
        return Monomial(exp_diff, self.cost_vector)

    def __eq__(self, other):
        return np.all(self.exponent == other.exponent)

    def __gt__(self, other):
        if (self.cost > other.cost):
            return True
        else:
            return self.grevlex_gt(other)

    def __repr__(self):
        return '[' + ','.join(map(str, self.exponent)) + ']'

    cdef bint grevlex_gt(self, Monomial other):
        cdef int i
        if self.degree > other.degree:
            return True
        elif other.degree > self.degree:
            return False
        else:
            for i in range(self.size-1, -1, -1):
                if other.exponent[i] > self.exponent[i]:
                    return True
                elif self.exponent[i] > other.exponent[i]:
                    return False
        return False

    cdef Monomial reduce_positive(self, Binomial other):
        pass

    cdef bint is_divisible_monomial(self, Monomial other):
        cdef int i
        for i in range(self.size):
            if self[i] < other[i]:
                return False
        return True

    cdef bint is_divisible_binomial(self, Binomial other):
        cdef int i
        for i in range(self.size):
            if other[i] > 0 and self[i] < other[i]:
                return False
        return True

cdef Monomial gcd(Monomial m1, Monomial m2):
    return Monomial(np.minimum(m1.exponent, m2.exponent), m1.cost_vector)

cdef Monomial lcm(Monomial m1, Monomial m2):
    return Monomial(np.maximum(m1.exponent, m2.exponent), m1.cost_vector)

"""
class Term {
public:
  friend Term operator*(const Term& t1, const Term& t2) { return Term{t1.coeff * t2.coeff, t1.monom * t2.monom}; }
  friend Term operator/(const Term& t1, const Term& t2) { return Term{t1.coeff / t2.coeff, t1.monom / t2.monom}; }
  friend bool operator==(const Term& t1, const Term& t2) { return (t1.coeff == t2.coeff) && (t1.monom == t2.monom); }
  friend bool operator!=(const Term& t1, const Term& t2) { return !(t1 == t2); }
  friend std::ostream& operator<<(std::ostream& os, const Term& t);

  Coefficient coeff;
  Monomial monom;
};
*/

class Binomial {
public:
  Binomial() : exponent{}, sug{0} {}
	Binomial(int dim) : exponent(dim,0), degree(0), cost(0) {}
  Binomial(std::initializer_list<int> exp, std::initializer_list<double> cost_vector);
  Binomial(std::vector<int> exp, std::vector<double> cost_vector);
  // Coefficient LC() const { return 1; }
  Monomial LM() const;
  int sugar() const { return sug; };
	bool is_nonzero() const;
  void flip();

  int operator[](int i) const { return exponent[i]; }
  int& operator[](int i) { return exponent[i]; }
	int size() const { return exponent.size();}

	double get_cost() const { return cost; }

	int get_degree() const { return degree; }
	void set_degree(int deg) {degree = deg;}
	
  //friend Polynomial operator+(const Polynomial& f1, const Polynomial& f2);
  //friend Polynomial operator-(const Polynomial& f1, const Polynomial& f2);
  //friend Polynomial operator*(const Term& t, const Polynomial& f);
  //friend Polynomial operator*(const Polynomial& f1, const Polynomial& f2);
	
	// Construct a new binomial given by difference of exponent vectors
	friend Binomial reduce_positive(const Binomial& f1, const Binomial& f2);
	// Construct a new binomial given by sum of exponent vectors
	friend Binomial reduce_negative(const Binomial& f1, const Binomial& f2);
	friend Binomial spoly(const Binomial& f, const Binomial& g);
	
	// Determine if positive part of f1 is larger than positive part of f2, if f1 divisible by f2
  friend bool is_divisible(const Binomial& f1, const Binomial& f2);
	// Determine if negative part of f1 is larger than positive part of f2, if f1 divisible by -f2
  friend bool is_divisible_negative(const Binomial& f1, const Binomial& f2); // Determine if 
  friend bool operator==(const Binomial& f1, const Binomial& f2);
  friend bool operator!=(const Binomial& f1, const Binomial& f2) { return !(f1 == f2); }
  friend std::ostream& operator<<(std::ostream& os, const Binomial& f);
private:
  std::vector<int> exponent; // Positive terms represent first monomial, negative the second monomial
  int sug;
  int degree;
  double cost;
};

/*
class Polynomial {
public:
  Polynomial() : terms{}, sug{0} {}
  Polynomial(std::initializer_list<Term> tms);
  Polynomial(std::vector<Term> tms);
  Coefficient LC() const { return terms[0].coeff; }
  Monomial LM() const { return terms[0].monom; }
  Term LT() const { return terms[0]; }
  int size() const { return terms.size(); };
  int sugar() const { return sug; };

  friend Polynomial operator+(const Polynomial& f1, const Polynomial& f2);
  friend Polynomial operator-(const Polynomial& f1, const Polynomial& f2);
  friend Polynomial operator*(const Term& t, const Polynomial& f);
  friend Polynomial operator*(const Polynomial& f1, const Polynomial& f2);
  friend bool operator==(const Polynomial& f1, const Polynomial& f2);
  friend bool operator!=(const Polynomial& f1, const Polynomial& f2) { return !(f1 == f2); }
  friend std::ostream& operator<<(std::ostream& os, const Polynomial& f);

  std::vector<Term> terms;

private:
  int sug;
};
*/

/**
 * Return the polynomial from the given string.
 *
 * Polynomials have the form
 *
 *  polynomial = term polynomial
 *        term = integer * monomial
 *             | integer
 *             | monomial
 *             | + term
 *             | - term
 *    monomial = variable ^ integer * monomial
 *             | variable ^ integer
 *             | variable * monomial
 *             | variable
 *
 * where variables are ["a", "b", ..., "h"], integers are sequences of
 * ["0", "1", ..., "9"], and no spaces are allowed.
 *
 * This function does not do any error checking, so be careful!
 */
//Polynomial parse_polynomial(const std::string& poly_string);


#endif


Monomial::Monomial(std::initializer_list<int> exp, std::initializer_list<double> cost_vector) {
  // if (exp.size() > N) throw std::invalid_argument("exponent vector is too large");
  std::fill(exponent.begin(), exponent.end(), 0);
  std::copy(exp.begin(), exp.end(), exponent.begin());
  degree = std::accumulate(exponent.begin(), exponent.end(), 0);
  cost = std::inner_product(exponent.begin(), exponent.end(), cost_vector.begin(), 0);
}


Monomial::Monomial(std::vector<int> exp, std::vector<double> cost_vector) {
	exponent = exp;
  degree = std::accumulate(exponent.begin(), exponent.end(), 0);
  cost = std::inner_product(exponent.begin(), exponent.end(), cost_vector.begin(), 0);
}  


Monomial operator*(const Monomial& m1, const Monomial& m2) {
  Monomial m(m1.size());
  m.degree = m1.degree + m2.degree;
  for (int i = 0; i < m1.size(); i++)
    m[i] = m1[i] + m2[i];
  return m;
}


Monomial operator/(const Monomial& m1, const Monomial& m2) {
  Monomial m(m1.size());
  m.degree = m1.degree - m2.degree;
  for (int i = 0; i < m1.size(); i++) {
    m[i] = m1[i] - m2[i];
  }
  return m;
}


bool operator>(const Monomial& m1, const Monomial& m2) {
	if (m1.cost > m2.cost)
		return true;
	else
		return grevlex_gt(m1, m2);
}

// grevlex
bool grevlex_gt(const Monomial& m1, const Monomial& m2) {
  if (m1.degree > m2.degree) {
    return true;
  } else if (m2.degree > m1.degree) {
    return false;
  } else {
    for (int i = m1.size()-1; i >= 0; i--) {
      if (m2.exponent[i] > m1.exponent[i])
	return true;
      else if (m1.exponent[i] > m2.exponent[i])
	return false;
    }
  }
  return false; 
}


bool operator==(const Monomial& m1, const Monomial& m2) {
  for (int i = 0; i < m1.size(); i++)
    if (m1.exponent[i] != m2.exponent[i]) return false;
  return true;
}


std::ostream& operator<<(std::ostream& os, const Monomial& m) {
  os << "[";
  for (int i = 0; i < m.size()-1; i++)
    os << m[i] << " ";
  os << m[m.size()-1] << "]";
  return os;
}

Monomial reduce_positive(const Monomial& f1, const Binomial& f2) {
	Monomial m(f1.size());
	for (int i = 0; i < f1.size(); i++)
		m[i] = f1[i] - f2[i];

	m.cost = f1.cost - f2.get_cost();
	m.degree = f1.degree - f2.get_degree();
	return m;
}

bool is_divisible(const Monomial& m1, const Monomial& m2) {
  for (int i = 0; i < m1.size(); i++) {
    if (m1[i] < m2[i]) return false;
  }
  return true;
}

bool is_divisible(const Monomial& m, const Binomial& b) {
  for (int i = 0; i < m.size(); i++) {
    if (b[i] > 0 && m[i] < b[i]) return false;
  }
  return true;
}


Monomial gcd(const Monomial& m1, const Monomial& m2) {
  Monomial m(m1.size());
  for (int i = 0; i < m1.size(); i++) {
    m[i] = std::min(m1[i], m2[i]);
    m.degree += m[i];
  }
  return m;
}


Monomial lcm(const Monomial& m1, const Monomial& m2) {
  Monomial m(m1.size());
  for (int i = 0; i < m1.size(); i++) {
    m[i] = std::max(m1[i], m2[i]);
    m.degree += m[i];
  }
  return m;
}

/*
std::ostream& operator<<(std::ostream& os, const Term& t) {
  if (t.monom == Monomial{}) {
    os << t.coeff;
    return os;
  }
  os << t.coeff << t.monom;
  return os;
}
*/

Monomial Binomial::LM() const{
	Monomial m(size());
	int degree = 0;
	for (int i = 0; i < size(); i++)
		if (exponent[i] > 0)
		{
			m[i] = exponent[i];
			degree += exponent[i];
		}
	m.set_degree(degree);
	return m;
}

Binomial::Binomial(std::initializer_list<int> exp, std::initializer_list<double> cost_vector) {
  // if (exp.size() > size()) throw std::invalid_argument("exponent vector is too large");
  std::fill(exponent.begin(), exponent.end(), 0);
  std::copy(exp.begin(), exp.end(), exponent.begin());
  degree = std::accumulate(exponent.begin(), exponent.end(), 0);
  cost = std::inner_product(exponent.begin(), exponent.end(), cost_vector.begin(), 0.0);
	if (cost < 0)
		flip();

  sug = LM().deg();
}

Binomial::Binomial(std::vector<int> exp, std::vector<double> cost_vector) {
	exponent = exp;
  degree = std::accumulate(exponent.begin(), exponent.end(), 0);
  cost = std::inner_product(exponent.begin(), exponent.end(), cost_vector.begin(), 0.0);
	if (cost < 0)
		flip();

  sug = LM().deg();
}  

bool Binomial::is_nonzero() const
{
	for (int i = 0; i < size(); i++)
		if (exponent[i] != 0)
			return true;
	return false;
}

Binomial reduce_positive(const Binomial& f1, const Binomial& f2) {
	Binomial b(f1.size());
	for (int i = 0; i < f1.size(); i++)
		b[i] = f1[i] - f2[i];

	b.degree = f1.degree - f2.degree;
	b.cost = f1.cost - f2.cost;
	return b;
}

Binomial reduce_negative(const Binomial& f1, const Binomial& f2) {
	Binomial b(f1.size());
	for (int i = 0; i < f1.size(); i++)
		b[i] = f1[i] + f2[i];

	b.degree = f1.degree + f2.degree;
	b.cost = f1.cost + f2.cost;
	return b;

}

bool is_divisible(const Binomial& f1, const Binomial& f2) {
  for (int i = 0; i < f1.size(); i++) {
    if (f2[i] > 0 && f2[i] > f1[i])
			return false;
  }
  return true;
}

bool is_divisible_negative(const Binomial& f1, const Binomial& f2) {
  for (int i = 0; i < f1.size(); i++) {
    if (f2[i] > 0 && f2[i] > -f1[i]) return false;
  }
  return true;

}

void Binomial::flip() {
	for (int i = 0; i < size(); i++)
		exponent[i] *= -1;
	cost *= -1;
}

/*
Polynomial operator+(const Polynomial& f1, const Polynomial& f2) {
  Polynomial g;
  g.sug = std::max(f1.sug, f2.sug);
  int i = 0, j = 0;
  while (i < f1.terms.size() && j < f2.terms.size()) {
    Term t1 = f1.terms[i];
    Term t2 = f2.terms[j];
    if (t1.monom > t2.monom) {
      g.terms.push_back(t1);
      i++;
    } else if (t2.monom > t1.monom) {
      g.terms.push_back(t2);
      j++;
    } else {
      Coefficient c = t1.coeff + t2.coeff;
      if (c != 0)
	g.terms.push_back({c, t1.monom});
      i++;
      j++;
    }
  }
  if (i < f1.terms.size()) {
    for (int k = i; k < f1.terms.size(); k++)
      g.terms.push_back(f1.terms[k]);
  } else {
    for (int k = j; k < f2.terms.size(); k++)
      g.terms.push_back(f2.terms[k]);
  }
  return g;
}
*/

/*
Polynomial operator-(const Polynomial& f1, const Polynomial& f2) {
  Polynomial f = f2;
  for (Term& t : f.terms)
    t.coeff = (-1) * t.coeff;
  return f1 + f;
}
*/

bool operator==(const Binomial& f1, const Binomial& f2) {
  for (int i = 0; i < f1.size(); i++)
    if (f1.exponent[i] != f2.exponent[i]) return false;
  return true;
}

/*
Polynomial operator*(const Term& t, const Polynomial& f) {
  Polynomial g;
  g.sug = t.monom.deg() + f.sug;
  for (const Term& ft : f.terms)
    g.terms.push_back(t * ft);
  return g;
}


Polynomial operator*(const Polynomial& f1, const Polynomial& f2) {
  Polynomial g;
  for (const Term& t : f1.terms)
    g = g + t * f2;
  return g;
}
*/

std::ostream& operator<<(std::ostream& os, const Binomial& f) {
  os << "[";
  for (int i = 0; i < f.size()-1; i++)
    os << f[i] << " ";
  os << f[f.size()-1] << "]";
  return os;
}

/*
int parse_variable(std::istringstream& is) {
  int i = is.get() - (int) 'a';
  if (i < 0 || i > N) throw std::invalid_argument("invalid variable name");
  return i;
}


Monomial parse_monomial(std::istringstream& is) {
  if (is.peek() == EOF) return Monomial{};
  int var = parse_variable(is);
  std::array<int, N> exp = {};
  switch(is.peek()) {
  case '^':
    is.get();
    int power;
    is >> power;
    exp[var] += power;
    if (is.peek() == '*') {
      is.get();
      return Monomial(exp) * parse_monomial(is);
    } else {
      return Monomial(exp);
    }
  case '*':
    is.get();
    exp[var] += 1;
    return Monomial(exp) * parse_monomial(is);
  default:
    exp[var] += 1;
    return Monomial(exp);
  }
}


Term parse_term(std::istringstream& is) {
  switch (is.peek()) {
  case '+':
    is.get();
    return parse_term(is);
  case '-':
    is.get();
    return Term{-1, {}} * parse_term(is);
  default:
    if (std::isdigit(is.peek())) {
      int c;
      is >> c;
      if (is.peek() == '*') {
	is.get();
	Monomial m = parse_monomial(is);
	return {Coefficient{c},  m};
      } else {
	return {Coefficient{c}, {}};
      }
    } else {
      Monomial m = parse_monomial(is);
      return {Coefficient{1},  m};
    }
  }
}


Polynomial parse_polynomial(std::istringstream& is) {
  if (is.peek() == EOF) {
    return Polynomial{};
  } else {
    Term t = parse_term(is);
    return Polynomial{t} + parse_polynomial(is);
  }
}


Polynomial parse_polynomial(const std::string& poly_string) {
  std::istringstream iss(poly_string);
  return parse_polynomial(iss);
}
*/

Binomial spoly(const Binomial& f, const Binomial& g) {
	Binomial s(f.size());
	for (int i = 0; i < f.size(); i++)
		s[i] = f[i] - g[i];

	s.cost = f.cost - g.cost;
	// All Binomials will by default correspond to improving moves
	if (s.cost < 0)
		s.flip();

  return s;
}
"""
