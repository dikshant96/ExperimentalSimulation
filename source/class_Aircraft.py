import numpy as np
class Wing:
    def __init__(self):
        self.b = 1.397
        self.S = 0.2172
        self.cmac = 0.165
        self.A = 8.98
        self.taper = 0.4
        self.dihedral = 4
        self.twist = 2
        self.tc = 0.15
        self.V = 0.003029

class Tailh:
    def __init__(self):
        self.b = 0.576
        self.S = 0.0858
        self.cmac = 0.149
        self.A = 3.87
        self.taper = 1
        self.dihedral = 6
        self.twist = 0
        self.tc = 0.15
        self.V = 0.0024485

class Tailv:
    def __init__(self):
        self.b = 0.258
        self.S = 0.0415
        self.cmac = 0.170
        self.A = 1.606
        self.taper = 0.43
        self.tc = 0.15
        self.V = 0.0003546

class Prop:
    def __init__(self):
        self.D = 0.2032
        self.Nblades = 6
        self.pitch = 45
        self.Sp = np.pi * (self.D/2)**2

class Fuselage:
    def __init__(self):
        self.l = 1.342
        self.d = 0.140
        self.V = 0.0160632


class Aircraft:
    def __init__(self):
        self.wing = Wing()
        self.tailh = Tailh()
        self.tailv = Tailv()
        self.prop = Prop()
        self.fuselage = Fuselage()