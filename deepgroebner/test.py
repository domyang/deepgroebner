import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

import geometric_buchberger as gb

cost = np.zeros(27)
cost_dist = gb.ConstantDistribution(cost)
x0 = np.random.randint(0, 6, 27)
x0_dist = gb.ConstantDistribution(x0)

markov_basis = np.loadtxt('333.1.mar', skiprows=1, dtype=int)
buchberger_env = gb.BuchbergerEnv(markov_basis, cost_dist, x0_dist, elimination='lcm')
buchberger_env.reset()

while buchberger_env.P.shape[0] > 0:
    print(buchberger_env.P[0])
    buchberger_env.step(0)
