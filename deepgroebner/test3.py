import argparse
import datetime
import gymnasium as gym
import json
import numpy as np
import networkx as nx
import os
import tensorflow as tf

from deepgroebner.buchberger import LeadMonomialsEnv, BuchbergerAgent
from deepgroebner.pg import PGAgent, PPOAgent
from deepgroebner.networks import MultilayerPerceptron, ParallelMultilayerPerceptron, AttentionPMLP, TransformerPMLP, PairsLeftBaseline, AgentBaseline, PointerNetwork, RecurrentValueModel, PoolingValueModel
from wrapped import CLeadMonomialsEnv
from deepgroebner.environments import VectorEnv, AlphabeticalEnv


import basis_methods.benchmark_problems as bp
import basis_methods.basis_construction as bc


