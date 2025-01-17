#include "polynomials.h"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <tuple>

/*
Coefficient operator/(Coefficient c1, Coefficient c2) {

  // compute a = c2 inverse using extended Euclidean algorithm
  int a = 0, a_ = 1;
  int b = P, b_ = c2.c;
  while (b_ != 0) {
    int q = b / b_;
    std::tie(a, a_) = std::make_tuple(a_, a - q * a_);
    std::tie(b, b_) = std::make_tuple(b_, b - q * b_);
  }

  return c1.c * a;
}
*/

Monomial::Monomial(std::initializer_list<int> exp, std::initializer_list<double> cost_vector) {
  // if (exp.size() > N) throw std::invalid_argument("exponent vector is too large");
  std::fill(exponent.begin(), exponent.end(), 0);
  std::copy(exp.begin(), exp.end(), exponent.begin());
  degree = std::accumulate(exponent.begin(), exponent.end(), 0);
  cost = std::inner_product(exponent.begin(), exponent.end(), cost_vector.begin(), 0.0);
}


Monomial::Monomial(std::vector<int> exp, std::vector<double> cost_vector) {
	exponent = exp;
  degree = std::accumulate(exponent.begin(), exponent.end(), 0);
  cost = std::inner_product(exponent.begin(), exponent.end(), cost_vector.begin(), 0.0);
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
  //std::cout << "exp size " << exp.size() << ", cost size " << cost_vector.size() << std::endl;
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
