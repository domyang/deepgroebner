import argparse
import datetime
import gymnasium as gym
import json
import numpy as np
import networkx as nx
import os
import tensorflow as tf

from pg import PGAgent, PPOAgent
from networks import MultilayerPerceptron, ParallelMultilayerPerceptron, AttentionPMLP, TransformerPMLP, PairsLeftBaseline, AgentBaseline, PointerNetwork, RecurrentValueModel, PoolingValueModel
from wrapped import CLeadMonomialsEnv


import basis_methods.benchmark_problems as bp
import basis_methods.basis_construction as bc


eps = float(0.2)
policy_lr = float(1e-4)
policy_updates = int(40)
value_lr = float(1e-3)
value_updates = int(40)
gam = float(0.99)
lam = float(0.97)
policy_kld_limit = float(0.01)
ent_bonus = float(0.0)
k = 2


n = 4
grid_nxn = nx.convert_node_labels_to_integers(nx.grid_2d_graph(n,n))
tv_mat_nxn = bp.tv_matrix_2d(n, with_bounds=True)
markov_basis = bc.markov_basis(tv_mat_nxn)


def tv_cost_vec(n_nodes, n_edges, c):
    return np.hstack([c, np.ones(n_edges*2), np.zeros(n_nodes)]).astype(np.float64)

n_nodes = len(grid_nxn.nodes)
n_edges = len(grid_nxn.edges)

def x0_vec(n_nodes, n_edges):
    return np.hstack([np.zeros(n_nodes + 2*n_edges), np.ones(n_nodes)]).astype(np.int32)

cost_vec = tv_cost_vec(n_nodes, n_edges, 2*np.random.randn(n_nodes))
cost_vec = cost_vec + 1e-8 * np.random.randn(cost_vec.size)
np.save('cost_vec.npy', cost_vec)
x0 = x0_vec(n_nodes, n_edges)


from distributions import Distribution

class TV_Cost_Dist(Distribution):
    def __init__(self, graph, cost, eps):
        self.graph = graph
        self.cost = cost
        self.eps = eps

    def sample(self):
        noise = np.zeros(2*len(self.graph.nodes) + 2*len(self.graph.edges))
        noise[:len(self.graph.nodes)] = self.eps*np.random.randn(len(self.graph.nodes))
        return self.cost + noise

tv_cost_dist = TV_Cost_Dist(grid_nxn, cost_vec, 0.5)

import pyomo.environ as pyenv

def tv_model(graph, cost):
    model = pyenv.ConcreteModel()
    model.nodes = list(graph.nodes)
    model.edges = list(graph.edges)
    model.c = pyenv.Param(model.nodes, initialize=cost, mutable=True)

    model.x = pyenv.Var(model.nodes, within=pyenv.NonNegativeReals, bounds=(0,1))
    model.a_plus = pyenv.Var(model.edges, within=pyenv.NonNegativeReals)
    model.a_minus = pyenv.Var(model.edges, within=pyenv.NonNegativeReals)

    model.constraints = pyenv.ConstraintList()
    for u, v in graph.edges:
        model.constraints.add(model.x[u] - model.x[v] == model.a_plus[u,v] - model.a_minus[u,v])

    model.obj = pyenv.Objective(expr=(sum(model.c[i]*model.x[i] for i in model.nodes)
                                      + sum(model.a_plus[u,v] + model.a_minus[u,v] for (u,v) in model.edges)),
                                sense=pyenv.minimize)
    return model

model = tv_model(grid_nxn, cost_vec[:len(grid_nxn.nodes)])

solver = pyenv.SolverFactory('gurobi')

def resolve(cost_vector):
    for i in model.nodes:
        model.c[i] = cost_vector[i]
    solver.solve(model)
    return pyenv.value(model.obj)


env = CLeadMonomialsEnv(markov_basis.astype(np.int32), x0, cost_vector=tv_cost_dist,
                        stop="Ops", max_ops=500, solver=resolve, base_vars=25)


state = env.reset()

episodes = 10
max_episode_length = 500
epochs = 1000
batch_size = 64
save_freq = 100


def make_logdir():
    """Return the directory name for this run."""
    run_name = 'run'
    time_string = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    run_name = time_string + '_' + run_name
    logdir = os.path.join('data/train', run_name)
    os.makedirs(logdir)
    return logdir


policy_network = ParallelMultilayerPerceptron(hidden_layers=[int(128)])
value_network = 'degree'

dim = env.basis().shape[1]

agent = PPOAgent(policy_network=policy_network, method='clip', eps=eps,
                 policy_lr=policy_lr, policy_updates=policy_updates,
                 value_network=None, value_lr=value_lr, value_updates=value_updates,
                 gam=gam, lam=lam, kld_limit=policy_kld_limit, ent_bonus=ent_bonus)
batch = np.zeros((1, 10, 3 * dim), dtype=np.float32)
policy_network(batch)

logdir = make_logdir()

agent.train(env, episodes=episodes, epochs=epochs, save_freq=save_freq, logdir=logdir,
            max_episode_length=max_episode_length, batch_size=batch_size)

from guppy import hpy
hp = hpy()
import pdb; pdb.set_trace()

