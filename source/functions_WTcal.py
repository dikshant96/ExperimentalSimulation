import numpy as np

class Oper:
    pass

def get_qtunnelcal(dPb,testsec):
    if testsec == 5:
        dPblim = [100,300]
        facQpB = [[0.51549, 2.32312, 4.62743e-05, 0.0],
                  [0.51549, 2.32312, 4.62743e-05, 0.0],
                  [0.51549, 2.32312, 4.62743e-05, 0.0]]
    elif testsec == 7:
        dPblim = [50,400]
        facQpB =[[0.909286, 2.44320, -3.79098e-05, 0.0],
                 [0.402538, 2.45530, -5.74383e-05, 0.0],
                 [3.101380, 2.42426, 5.22943e-05, 0.0]]
    elif testsec == 9:
        dPblim = [100,5000]
        facQpB = [[0.0193846, 2.33053, 1.72980e-4, 0.0],
                  [0.1903500, 2.33534, 5.07273e-5, 0.0],
                  [0.1903500, 2.33534, 5.07273e-5, 0.0]]
    qInf = np.zeros((len(dPb),))
    for i,q in enumerate(qInf):
        if dPb[i] < dPblim[0]:
            idxQ = 0
        elif dPb[i] < dPblim[1]:
            idxQ = 1
        else:
            idxQ = 2
        qInf[i] = facQpB[idxQ][0]*dPb[i]**0 + facQpB[idxQ][1]*dPb[i]**1 + facQpB[idxQ][2]*dPb[i]**2 + facQpB[idxQ][3]*dPb[i]**3
    return qInf

def get_oper(BAL_data, testsec):
    oper = Oper()
    data = BAL_data.data
    index = BAL_data.index
    oper.dPb = data[:,index.dPb]
    oper.tInf = data[:,index.temp]+273.15
    oper.pBar = data[:,index.pBar]*100
    oper.AoA = data[:,index.AoA]

    oper.qInf = get_qtunnelcal(oper.dPb, testsec)
    oper.pInf = oper.pBar
    oper.rho = oper.pInf/(oper.tInf*287.05)
    oper.vInf = np.sqrt(2*oper.qInf/oper.rho)

    oper.mu = 1.716e-5*(oper.tInf/273.15)**1.5*(273.15+110.4)/(oper.tInf+110.4)
    oper.nu = oper.mu/oper.rho

    return oper