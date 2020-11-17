#ifndef BUCHBERGER_H
#define BUCHBERGER_H

#include "polynomials.h"

#include <map>
#include <queue>
#include <vector>


Polynomial spoly(const Polynomial& f, const Polynomial& g);
Polynomial reduce(const Polynomial& g, const std::vector<Polynomial>& F);
std::vector<Polynomial> minimalize(std::vector<Polynomial>& G);
std::vector<Polynomial> interreduce(std::vector<Polynomial>& G);
std::vector<Polynomial> buchberger(const std::vector<Polynomial>& F);

struct SPair {
  int i;
  int j;
  bool valid = true;
};

class SPairSet {
public:
  SPairSet(std::vector<Polynomial>& gens) : F{gens} {}
  int size() { return sz; }
  bool empty() const { return sz == 0; }
  SPair pop();
  void update(const Polynomial& r);

private:
  int sz = 0;
  std::vector<Polynomial>& F;
  std::priority_queue<std::tuple<int, int, int>> keys;
  std::map<std::tuple<int, int, int>, std::vector<SPair>> bins;
};

class BuchbergerEnv {
public:
  BuchbergerEnv() : state{0} {}
  BuchbergerEnv(int s) : state{s} {}
  int reset() { return state; }
  int step(int i, int j) { return i + j; }
private:
  int state;
};

#endif
