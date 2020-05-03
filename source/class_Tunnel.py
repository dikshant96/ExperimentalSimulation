import numpy as np

class Tunnel:
    def __init__(self):
        self.b = 1800 / 1000.
        self.h = 1250 / 1000.
        self.C = self.b * self.h - (0.5 * 300 * 300 * 4) / 1000000.
        self.L = self.h/self.b
