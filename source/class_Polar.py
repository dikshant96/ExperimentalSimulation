import numpy as np
from math import radians as rad

class Point:
    def __init__(self,alpha,beta,U,rho,qInf,temp,pInf,nu,n,F1,F2,F3,M1,M2,M3,aircraft):
        self.alpha = alpha
        self.beta = beta
        self.U = U
        self.rho = rho
        self.qInf = qInf
        self.temp = temp
        self.pInf = pInf
        self.nu = nu
        self.n = n
        self.F1 = F1
        self.F2 = F2
        self.F3 = F3
        self.M1 = M1
        self.M2 = M2
        self.M3 = M3
        self.get_reynolds(aircraft)
        self.get_advanceratio(aircraft)
        self.get_coeffs(aircraft)
        self.get_modelaxiscoeffs(aircraft)
        self.get_conventionalcoeffs(aircraft)

    def get_reynolds(self,aircraft):
        self.Re = (self.U*aircraft.wing.cmac)/self.nu

    def get_advanceratio(self,aircraft):
        if np.abs(self.n) < 1e-2:
            self.J = 0
        else:
            self.J = self.U/(self.n*aircraft.prop.D)

    def get_coeffs(self,aircraft):
        self.CF1 = self.F1 / (self.qInf * aircraft.wing.S)
        self.CF2 = self.F2 / (self.qInf * aircraft.wing.S)
        self.CF3 = self.F3 / (self.qInf * aircraft.wing.S)
        self.CM1 = self.M1 / (self.qInf * aircraft.wing.S * aircraft.wing.b)
        self.CM2 = self.M2 / (self.qInf * aircraft.wing.S * aircraft.wing.cmac)
        self.CM3 = self.M3 / (self.qInf * aircraft.wing.S * aircraft.wing.b)

    def get_modelaxiscoeffs(self,aircraft):
        self.CFt = self.CF1 * np.cos(rad(self.alpha)) - self.CF3 * np.sin(rad(self.alpha))
        self.CFn = self.CF3 * np.cos(rad(self.alpha)) + self.CF1 * np.sin(rad(self.alpha))
        self.CFs = self.CF2

        self.CMr = self.CM1 + (0.0465 / aircraft.wing.cmac)*self.CF2
        self.CMp = self.CM2 - (0.0465 / aircraft.wing.cmac)*self.CF1
        self.CMy = self.CM3

        self.CMr = self.CMr * np.cos(rad(self.alpha)) - self.CMy * np.sin(rad(self.alpha))
        self.CMp = self.CMp
        self.CMy = self.CMy * np.cos(rad(self.alpha)) + self.CMr * np.sin(rad(self.alpha))

    def get_conventionalcoeffs(self,aircraft):
        self.CFl = self.CFn * np.cos(rad(self.alpha)) - self.CFt * np.sin(rad(self.alpha))
        self.CFd = self.CFn * np.sin(rad(self.alpha)) + self.CFt * np.cos(rad(self.alpha))
        self.CFy = self.CFs
        self.CMp25c = self.CMp

class Polar:
    def __init__(self,raw_data,raw_index,aircraft):
        points = []
        for i in range(len(raw_data[:,0])):
            raw_alpha = raw_data[i,raw_index.AoA]
            raw_beta = raw_data[i, raw_index.AoS]
            raw_U = raw_data[i, raw_index.V]
            raw_rho = raw_data[i, raw_index.rho]
            raw_qInf = raw_data[i, raw_index.q]
            raw_temp = raw_data[i, raw_index.temp]
            raw_pInf = raw_data[i, raw_index.pInf]
            raw_nu = raw_data[i, raw_index.nu]
            raw_n = raw_data[i, raw_index.rpsM2]
            raw_F1 = raw_data[i, raw_index.F1]
            raw_F2 = raw_data[i, raw_index.F2]
            raw_F3 = raw_data[i, raw_index.F3]
            raw_M1 = raw_data[i, raw_index.M1]
            raw_M2 = raw_data[i, raw_index.M2]
            raw_M3 = raw_data[i, raw_index.M3]
            points.append(Point(raw_alpha,raw_beta,raw_U,raw_rho,raw_qInf,raw_temp,
                                raw_pInf,raw_nu,raw_n,raw_F1,raw_F2,raw_F3,raw_M1,
                                raw_M2,raw_M3,aircraft))
        self.points = points

