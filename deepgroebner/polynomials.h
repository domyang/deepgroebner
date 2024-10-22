#ifndef POLYNOMIALS_H
#define POLYNOMIALS_H

#include <array>
#include <iostream>
#include <string>
#include <vector>

/*
constexpr int P = 32003;
class Coefficient {
public:
  Coefficient() : c{0} {}
  Coefficient(int i) : c{(i < 0) ? (i % P) + P : i % P} {}

  friend Coefficient operator+(Coefficient c1, Coefficient c2) { return c1.c + c2.c; }
  friend Coefficient operator-(Coefficient c1, Coefficient c2) { return c1.c - c2.c; }
  friend Coefficient operator*(Coefficient c1, Coefficient c2) { return c1.c * c2.c; }
  friend Coefficient operator/(Coefficient c1, Coefficient c2);
  friend bool operator==(Coefficient c1, Coefficient c2) { return c1.c == c2.c; }
  friend bool operator!=(Coefficient c1, Coefficient c2) { return c1.c != c2.c; }
  friend std::ostream& operator<<(std::ostream& os, Coefficient c) { return os << c.c; }

private:
  int c;
};
*/
class Binomial;

class Monomial {
public:
  Monomial() : exponent{}, degree{0}, cost{0} {}
	Monomial(int dim) : exponent(dim,0), degree(0), cost(0) {}
  Monomial(std::initializer_list<int> exp, std::initializer_list<double> cost_vector);
  Monomial(std::vector<int> exp) : exponent{exp}, degree{0}, cost{0} {}
  Monomial(std::vector<int> exp, std::vector<double> cost_vector);

  int operator[](int i) const { return exponent[i]; }
  int& operator[](int i) { return exponent[i]; }
  int deg() const { return degree; }
	int size() const { return exponent.size();}
	double get_cost() const { return cost; }

	int get_degree() const { return degree; }
	void set_degree(int deg) {degree = deg;}

  friend Monomial operator*(const Monomial& m1, const Monomial& m2); // addition of degree vectors
  friend Monomial operator/(const Monomial& m1, const Monomial& m2); // subtraction of degree vectors
	friend Monomial reduce_positive(const Monomial& f1, const Binomial& f2);
  friend bool operator==(const Monomial& m1, const Monomial& m2);
  friend bool operator!=(const Monomial& m1, const Monomial& m2) { return !(m1 == m2); }
  friend bool operator>(const Monomial& m1, const Monomial& m2);
  friend bool grevlex_gt(const Monomial& m1, const Monomial& m2);
  friend bool operator<(const Monomial& m1, const Monomial& m2) { return m2 > m1; }
  friend std::ostream& operator<<(std::ostream& os, const Monomial& m);

  friend bool is_divisible(const Monomial& m1, const Monomial& m2);
  friend bool is_divisible(const Monomial& m1, const Binomial& m2);
  friend Monomial gcd(const Monomial& m1, const Monomial& m2); // elementwise minimum of degree
  friend Monomial lcm(const Monomial& m1, const Monomial& m2); // elementwise maximum of degree

private:
  std::vector<int> exponent;
  int degree;
  double cost;
};

/*
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
