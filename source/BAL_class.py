import numpy as np
import source.functions_WTcal as fWT
import source.functions_BALcal as fB

class Index:
    def __init__(self):
        self.run = 0
        self.hr = 1
        self.min = 2
        self.sec = 3
        self.AoA = 4
        self.AoS = 5
        self.dPb = 6
        self.pBar = 7
        self.temp = 8
        self.B1 = 9
        self.B2 = 10
        self.B3 = 11
        self.B4 = 12
        self.B5 = 13
        self.B6 = 14
        self.rpmWT = 15
        self.rho = 16
        self.q = 17
        self.V = 18
        self.Re = 19
        self.rpsM1 = 20
        self.rpsM2 = 21
        self.iM1 = 22
        self.iM2 = 23
        self.dPtQ = 24

        self.pInf = 25
        self.nu = 26
        self.B1_cal = 27
        self.B2_cal = 28
        self.B3_cal = 29
        self.B4_cal = 30
        self.B5_cal = 31
        self.B6_cal = 32
        self.F1 = 33
        self.F2 = 34
        self.F3 = 35
        self.M1 = 36
        self.M2 = 37
        self.M3 = 38


class Index_0:
    def __init__(self):
        self.run = 0
        self.hr = 1
        self.min = 2
        self.sec = 3
        self.AoA = 4
        self.AoS = 5
        self.pBar = 6
        self.temp = 7
        self.B1 = 8
        self.B2 = 9
        self.B3 = 10
        self.B4 = 11
        self.B5 = 12
        self.B6 = 13


class BAL_data:
    def __init__(self,fname,fname0,testsec):
        self.fname = fname
        self.index = Index()
        self.index0 = Index_0()
        self.data = self.read_data()
        self.fname0 = fname0
        self.data0 = self.read_zerodata()

        'Process WT data'
        self.data = np.append(self.data,np.zeros((len(self.data),2)), axis=1)
        self.process_WTcal(testsec)

        'Process cal data'
        data_cal = fB.process_BALcallibration(self)
        self.data = np.append(self.data,data_cal,axis=1)


    def read_data(self):
        f = open(self.fname)
        lines = f.readlines()
        f.close()
        data = []
        dPbCut = 25.0
        for i, line in enumerate(lines):
            if i > 1:
                line = line.strip()
                line = line.replace(':', '\t')
                line = line.split()
                line = line[:25]
                line = list(map(lambda x:float(x), line))
                if line[self.index.dPb] > dPbCut:
                   data.append(line)
        data = np.array(data)
        return data

    def read_zerodata(self):
        f = open(self.fname0)
        lines = f.readlines()
        f.close()
        data = []
        for i, line in enumerate(lines):
            if i > 1:
                line = line.strip()
                line = line.replace(':', '\t')
                line = line.split()
                line = line[:14]
                line = list(map(lambda x: float(x), line))
                data.append(line)

        data = np.array(data)
        return data

    def process_WTcal(self,testsec):
        oper = fWT.get_oper(self, testsec)
        for i in range(len(self.data)):
            self.data[i, self.index.q] = oper.qInf[i]
            self.data[i, self.index.rho] = oper.rho[i]
            self.data[i, self.index.V] = oper.vInf[i]
            self.data[i, self.index.temp] = oper.tInf[i]
            self.data[i, self.index.pInf] = oper.pInf[i]
            self.data[i, self.index.pBar] = oper.pBar[i]
            self.data[i, self.index.nu] = oper.nu[i]
