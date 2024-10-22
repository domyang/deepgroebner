import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import graphviz as gv
import networkx as nx
import pyomo.environ as pyenv

import basis_methods.benchmark_problems as bp
import basis_methods.basis_construction as bc

from wrapped import CLeadMonomialsEnv
from distributions import Distribution

def vc_x0(graph):
    return np.hstack([np.ones(len(graph.nodes) + len(graph.edges)), np.zeros(len(graph.nodes))]).astype(np.int32)

grid_5x5 = nx.convert_node_labels_to_integers(nx.grid_2d_graph(5,5))

design_mat = bp.vc_matrix(grid_5x5, with_bounds=True)
markov_basis = bc.markov_basis(design_mat)
x0 = vc_x0(grid_5x5)

dim = markov_basis.shape[1]

print("Creating Environment")
cost_vector = np.zeros(dim, dtype=np.float64)
cost_vector[:len(grid_5x5)] = 1



def vertex_cover_problem(graph, costs):
    model = pyenv.ConcreteModel()
    model.graph = graph
    model.x = pyenv.Var(list(graph.nodes), within=pyenv.Binary)
    model.c = pyenv.Param(list(graph.nodes), initialize=costs, mutable=True)
    model.cover_constraints = pyenv.ConstraintList()
    for i, j in graph.edges:
        model.cover_constraints.add(model.x[i] + model.x[j] >= 1)
    model.obj = pyenv.Objective(expr=sum(model.c[i] * model.x[i] for i in graph.nodes), sense=pyenv.minimize)
    return model

model = vertex_cover_problem(grid_5x5, cost_vector[:25].astype(int))

solver = pyenv.SolverFactory('gurobi')

def resolve(costs):
    for i in model.graph.nodes:
        model.c[i] = costs[i]

    solver.solve(model)
    return pyenv.value(model.obj)


class VCDistribution(Distribution):
    def __init__(self, graph, cost, eps):
        self.graph = graph
        self.cost = cost
        self.eps = eps
    def sample(self):
        noise = np.zeros(len(self.graph.nodes)*2 + len(self.graph.edges))
        noise[:len(self.graph.nodes)] = np.random.randn(len(self.graph.nodes))*self.eps
        return self.cost + noise

cost_dist = VCDistribution(grid_5x5, cost_vector, 0.5)

env = CLeadMonomialsEnv(markov_basis.astype(np.int32), x0, cost_dist, solver=resolve)
print(env.reset())

print(env.value())
