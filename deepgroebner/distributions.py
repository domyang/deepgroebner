import numpy as np

class Distribution:
    def sample(self):
        pass


class ConstantDistribution(Distribution):
    def __init__(self, point):
        self.point = point

    def sample(self):
        return self.point
