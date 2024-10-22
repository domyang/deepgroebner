# distutils: language = c++
# distutils: sources = deepgroebner/polynomials.cpp
# distutils: extra_compile_args = -std=c++17
# cython: language_level = 3

import numpy as np

from libcpp.vector cimport vector
from .buchberger cimport LeadMonomialsEnv, Numpy1DArray

from distributions import Distribution, ConstantDistribution

cdef class CLeadMonomialsEnv:

    cdef dict __dict__
    cdef LeadMonomialsEnv c_env

    def __cinit__(self, markov_basis, x0, cost_vector, stop="Ops", max_ops=500, solver=None, base_vars=None):
        cdef Numpy1DArray ar_descr
        cdef vector[Numpy1DArray] vec
        cdef int[:] ar

        dimension = markov_basis.shape[1]
        for ar in markov_basis:
            ar_descr.size = ar.size
            if ar.size > 0:
                ar_descr.ptr = &ar[0]
            else:
                ar_descr.ptr = NULL
            vec.push_back(ar_descr)

        if isinstance(x0, Distribution):
            x0_sample = x0.sample()
            self.x0_dist = x0
        else:
            x0_sample = x0
            self.x0_dist = ConstantDistribution(x0)

        if isinstance(cost_vector, Distribution):
            cost_sample = cost_vector.sample()
            self.cost_dist = cost_vector
        else:
            cost_sample = cost_vector
            self.cost_dist = ConstantDistribution(cost_vector)

        cdef int[:] x0_arr = x0_sample
        cdef double[:] cost_arr = cost_sample

        self.solver = solver
        obj = solver(cost_sample)
        self.obj = obj
        self.base_vars = base_vars
        self.x0 = x0
        self.cost = cost_sample

        self.c_env = LeadMonomialsEnv(vec, &x0_arr[0], False, True,
                                      dimension, &cost_arr[0], stop.encode(), max_ops, obj)

    def reset(self):
        x0_sample = self.x0_dist.sample()
        cost_sample = self.cost_dist.sample()
        self.cost
        cdef int[:] x0_arr = x0_sample
        cdef double[:] cost_arr = cost_sample

        obj = self.solver(cost_sample)
        self.obj = obj

        self.c_env.reset(&x0_arr[0], &cost_arr[0], obj)
        state = np.array(self.c_env.state, dtype=np.int32).reshape(-1, self.c_env.cols)
        #if self.base_vars:
        #    state1 = state[:,:self.base_vars]
        #    state2 = state[:,len(self.x0):len(self.x0)+self.base_vars]
        #    cost = self.cost[:self.base_vars]
        #    state = np.hstack([state1, state2, np.tile(cost, (state1.shape[0],1))])
        #else:
        state = np.hstack([state, np.tile(self.cost, (state.shape[0],1))])
        return state

    def step(self, action):
        #print(f"Step {action}, {self.pairs()[action]}")
        reward, stop = self.c_env.step(action)
        state = np.array(self.c_env.state, dtype=np.int32).reshape(-1, self.c_env.cols)
        #if self.base_vars:
        #    state1 = state[:,:self.base_vars]
        #    state2 = state[:,len(self.x0):len(self.x0)+self.base_vars]
        #    cost = self.cost[:self.base_vars]
        #    state = np.hstack([state1, state2, np.tile(cost, (state1.shape[0],1))])
        #else:
        state = np.hstack([state, np.tile(self.cost, (state.shape[0],1))])
        return state, reward, stop, {}

    def seed(self, seed=None):
        if seed is not None:
            self.c_env.seed(seed)

    def basis(self):
        basis = np.array(self.c_env.basis, dtype=np.int32).reshape(-1, self.c_env.cols // 2)
        return basis

    def pairs(self):
        pairs = np.array(self.c_env.pairs, dtype=np.int32).reshape(-1, 2)
        return pairs

    def select(self, strategy):
        return self.c_env.select(strategy)

    def complete(self, strategy=0):
        while len(self.pairs() > 0):
            idx = self.select(strategy)
            _, reward, stop, _ = self.step(idx)
            if stop:
                break

    def get_cost(self):
        return np.array(self.c_env.get_cost(), dtype=np.float64)

    def objective(self):
        return self.c_env.objective()
    
    def value(self, strategy='degree', gamma=0.99, max_iters=500):
        return self.c_env.value(strategy.encode(), gamma, max_iters)

    def copy(self):
        copy = CLeadMonomialsEnv()
        copy.c_env = LeadMonomialsEnv(self.c_env)
        return copy
