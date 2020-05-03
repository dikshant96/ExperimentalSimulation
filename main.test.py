import source.BAL_class
import pandas as pd
import numpy as np
import source.functions_WTcal
import source.class_Aircraft as cA
import source.class_Polar as cP
import source.class_Tunnel as cT
import matplotlib.pyplot as plt
import source.functions_blockage as fb
import copy

aircraft = cA.Aircraft()
tunnel = cT.Tunnel()

'Tail Off case'
dir = 'data/'
polar_tailOff = dir+'raw_tailOff_bal.txt'
polar_tailOff_0 = dir+'zero_tailOff_bal.txt'

polar_tailOff_raw = source.BAL_class.BAL_data(polar_tailOff,polar_tailOff_0,5)
raw_data_40 = []
for data_i in polar_tailOff_raw.data:
    if np.abs(data_i[polar_tailOff_raw.index.V]-40) < 2:
        raw_data_40.append(data_i[:])
raw_data_40 = np.array(raw_data_40)

polar_tailOff_uncorr = cP.Polar(raw_data_40,polar_tailOff_raw.index,aircraft)
polar_tailOff_corr = copy.deepcopy(polar_tailOff_uncorr)
fb.correct_notail(polar_tailOff_corr,aircraft,tunnel,'w')

CLu_tailOff_array = fb.get_listfrompolar(polar_tailOff_uncorr,'CL')
CDu_tailOff_array = fb.get_listfrompolar(polar_tailOff_uncorr,'CD')
CMu_tailOff_array = fb.get_listfrompolar(polar_tailOff_uncorr,'CM')
CLc_tailOff_array = fb.get_listfrompolar(polar_tailOff_corr,'CL')
CDc_tailOff_array = fb.get_listfrompolar(polar_tailOff_corr,'CD')
CMc_tailOff_array = fb.get_listfrompolar(polar_tailOff_corr,'CM')

'Tail On case'
polar_tailOn = dir+'raw_tailOn_bal.txt'
polar_tailOn_0 = dir+'zero_tailOn_bal.txt'
polar_tailOn_raw = source.BAL_class.BAL_data(polar_tailOn,polar_tailOn_0,5)
polar_tailOn_uncorr = cP.Polar(polar_tailOn_raw.data,polar_tailOn_raw.index,aircraft)
polar_tailOn_corr = copy.deepcopy(polar_tailOn_uncorr)
fb.correct_tail(polar_tailOn_corr,polar_tailOff_uncorr,aircraft,tunnel,'wt')

CLu_tailOn_array = fb.get_listfrompolar(polar_tailOn_uncorr,'CL')
CDu_tailOn_array = fb.get_listfrompolar(polar_tailOn_uncorr,'CD')
CMu_tailOn_array = fb.get_listfrompolar(polar_tailOn_uncorr,'CM')
CLc_tailOn_array = fb.get_listfrompolar(polar_tailOn_corr,'CL')
CDc_tailOn_array = fb.get_listfrompolar(polar_tailOn_corr,'CD')
CMc_tailOn_array = fb.get_listfrompolar(polar_tailOn_corr,'CM')
CLwc_tailOn_array = fb.get_listfrompolar(polar_tailOn_corr,'CLw')
CDwc_tailOn_array = fb.get_listfrompolar(polar_tailOn_corr,'CDw')
CMwc_tailOn_array = fb.get_listfrompolar(polar_tailOn_corr,'CMw')
CLtc_tailOn_array = fb.get_listfrompolar(polar_tailOn_corr,'CLt')
CDtc_tailOn_array = fb.get_listfrompolar(polar_tailOn_corr,'CDt')
CMtc_tailOn_array = fb.get_listfrompolar(polar_tailOn_corr,'CMt')

'Prop On Case'
polar_J16 = dir+'raw_J16_PropOn_bal.txt'
polar_J16_0 = dir+'zero_PropOn_bal.txt'
polar_J16_raw = source.BAL_class.BAL_data(polar_J16,polar_J16_0,5)
polar_J16_uncorr = cP.Polar(polar_J16_raw.data,polar_J16_raw.index,aircraft)
fb.get_thrust_curve(polar_J16_uncorr,aircraft)

alpha_tailOn =np.arange(-4,14+2,2)
Lwh = np.array([point.F3 for point in polar_tailOn_uncorr.points])
Dwh = np.array([point.F1 for point in polar_tailOn_uncorr.points])
#Mwh = np.array([point.M2 for point in polar_tailOn_uncorr.points])

alpha_tailOff = []
Lw = []
Dw = []
for i in range(len(alpha_tailOn)):
    alpha_u = np.array([point.alpha for point in polar_tailOff_uncorr.points])
    Lw_u = np.array([point.F3 for point in polar_tailOff_uncorr.points])
    Dw_u = np.array([point.F1 for point in polar_tailOff_uncorr.points])
    min_index = np.argmin(np.abs(alpha_u-alpha_tailOn[i]))
    alpha_tailOff.append(alpha_u[min_index])
    Lw.append(Lw_u[min_index])
    Dw.append(Dw_u[min_index])

alpha_tailOff = np.array(alpha_tailOff)
Lw = np.array(Lw)
Dw = np.array(Dw)
Lh = Lwh - Lw
Dh = Dwh - Dw

alpha_prop = np.arange(-4,8+4,4)
Fz = np.array([point.F3 for point in polar_J16_uncorr.points])
Fx = np.array([point.F1 for point in polar_J16_uncorr.points])
qInf = np.array([point.qInf for point in polar_J16_uncorr.points])

Lw_prop = Lw[0:8:2]
Dw_prop = Dw[0:8:2]
Lh_prop = Lh[0:8:2]
Dh_prop = Dh[0:8:2]
Lwh_prop = Lwh[0:8:2]
Dwh_prop = Dwh[0:8:2]


CLw_prop = Lw_prop/(qInf*aircraft.wing.S)
CDw_prop = Dw_prop/(qInf*aircraft.wing.S)
CLh_prop = Lh_prop/(qInf*aircraft.wing.S)
CDh_prop = Dh_prop/(qInf*aircraft.wing.S)
CLwh_prop = Lwh_prop/(qInf*aircraft.wing.S)
CDwh_prop = Dwh_prop/(qInf*aircraft.wing.S)

CLh_alpha, c = np.polyfit(alpha_prop,CLh_prop,1)
dCdh_dClh2, c2 = np.polyfit(CLh_prop**2,CDh_prop,1)

CLh_fit = CLh_alpha*alpha_prop + c
CDh_fit = dCdh_dClh2*CLh_prop**2 + c2
CFz = Fz/(qInf*aircraft.wing.S)
plt.plot(alpha_prop, CDh_prop)

Tc = polar_J16_uncorr.points[1].Tc
q_ratio = fb.get_slipstream(Tc,aircraft)
delta_CLh_slipstream = CLh_prop[1]*(q_ratio - 1)
print(delta_CLh_slipstream)
delta_CLh_net = CFz[1] - CLwh_prop[1]
delta_CLh_induction = delta_CLh_net - delta_CLh_slipstream
CLh_alpha_b = CLh_alpha
delta_alpha_induction = delta_CLh_induction/CLh_alpha_b
#delta_drag_slipstream =