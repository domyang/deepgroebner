"""Tests for policy gradient agents."""

import numpy as np
import pytest

from deepgroebner.pg import *


@pytest.mark.parametrize("rewards, gam, discounted_rewards", [
    ([], 0.9, []),
    ([1, 2, 3], 1, [6, 5, 3]),
    ([1, 1, 1, 1], 0.9, [3.439, 2.71, 1.9, 1.]),
    (np.array([]), 0.9, np.array([])),
    (np.array([1, 2, 3]), 1, np.array([6, 5, 3])),
    (np.array([1., 1., 1., 1.]), 0.9, np.array([3.439, 2.71, 1.9, 1.])),
])
def test_discount_rewards_0(rewards, gam, discounted_rewards):
    assert np.array_equal(discount_rewards(rewards, gam), discounted_rewards)


def test_discount_rewards_1():
    L = [1, 2, 3, 4, 5]
    L[2:] = discount_rewards(L[2:], 0.5)
    assert L == [1, 2, 6.25, 6.5, 5]


@pytest.mark.parametrize("gam, lam, advantages", [
    (1.0, 1.0, np.array([5., 4., 3., 2., 1.])),
    (0.5, 1.0, np.array([1.9375, 1.875, 1.75, 1.5, 1.])),
    (1.0, 0.5, np.array([1.9375, 1.875, 1.75, 1.5, 1.])),
    (0.5, 0.5, np.array([1.33203125, 1.328125, 1.3125, 1.25, 1.])),
])
def test_compute_advantages_0(gam, lam, advantages):
    rewards = [1, 1, 1, 1, 1]
    values = [0, 0, 0, 0, 0]
    assert np.array_equal(compute_advantages(rewards, values, gam, lam), advantages)
