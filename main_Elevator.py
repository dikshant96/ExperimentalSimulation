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

CLu_tailOff = fb.get_listfrompolar(polar_tailOff_uncorr,'CL')
CDu_tailOff = fb.get_listfrompolar(polar_tailOff_uncorr,'CD')
CMu_tailOff = fb.get_listfrompolar(polar_tailOff_uncorr,'CM')
CLc_tailOff = fb.get_listfrompolar(polar_tailOff_corr,'CL')
CDc_tailOff = fb.get_listfrompolar(polar_tailOff_corr,'CD')
CMc_tailOff = fb.get_listfrompolar(polar_tailOff_corr,'CM')

'Tail On case'
polar_tailOn = dir+'raw_tailOn_bal.txt'
polar_tailOn_0 = dir+'zero_tailOn_bal.txt'
polar_tailOn_raw = source.BAL_class.BAL_data(polar_tailOn,polar_tailOn_0,5)
polar_tailOn_uncorr = cP.Polar(polar_tailOn_raw.data,polar_tailOn_raw.index,aircraft)
polar_tailOn_corr = copy.deepcopy(polar_tailOn_uncorr)
fb.correct_tail(polar_tailOn_corr,polar_tailOff_uncorr,aircraft,tunnel,'wt')

CLu_tailOn = fb.get_listfrompolar(polar_tailOn_uncorr,'CL')
CDu_tailOn = fb.get_listfrompolar(polar_tailOn_uncorr,'CD')
CMu_tailOn = fb.get_listfrompolar(polar_tailOn_uncorr,'CM')
CLc_tailOn = fb.get_listfrompolar(polar_tailOn_corr,'CL')
CDc_tailOn = fb.get_listfrompolar(polar_tailOn_corr,'CD')
CMc_tailOn = fb.get_listfrompolar(polar_tailOn_corr,'CM')
CLwc_tailOn = fb.get_listfrompolar(polar_tailOn_corr,'CLw')
CDwc_tailOn = fb.get_listfrompolar(polar_tailOn_corr,'CDw')
CMwc_tailOn = fb.get_listfrompolar(polar_tailOn_corr,'CMw')
CLtc_tailOn = fb.get_listfrompolar(polar_tailOn_corr,'CLt')
CDtc_tailOn = fb.get_listfrompolar(polar_tailOn_corr,'CDt')
CMtc_tailOn = fb.get_listfrompolar(polar_tailOn_corr,'CMt')

CDtc_tailOn[:,1] = CDtc_tailOn[:,1]*(aircraft.wing.S/aircraft.tailh.S)
CLtc_tailOn[:,1] = CLtc_tailOn[:,1]*(aircraft.wing.S/aircraft.tailh.S)
dCD_CL2_tail = fb.get_tail_CDCL2(CDtc_tailOn,CLtc_tailOn,aircraft)

J16_alpha = [0.0, 4.0]
J16_Tc = [0.318,0.324]
J16_Tc_fit = np.interp(3,J16_alpha,J16_Tc)

J20_alpha = [0.0, 4.0]
J20_Tc = [0.1015,0.1127]
J20_Tc_fit = np.interp(3,J20_alpha,J20_Tc)

J24_alpha = [0.0, 4.0]
J24_Tc = [-0.0117,-0.0056]
J24_Tc_fit = np.interp(3,J24_alpha,J24_Tc)

Juc = np.array([1.6,2.0,2.4])
Tc = np.array([J16_Tc_fit,J20_Tc_fit,J24_Tc_fit])

'Prop On Cases'
polar_J16 = dir+'raw_J16_PropOn_bal.txt'
polar_J16_0 = dir+'zero_PropOn_bal.txt'
polar_J16_raw = source.BAL_class.BAL_data(polar_J16,polar_J16_0,5)
polar_J16_uncorr = cP.Polar(polar_J16_raw.data,polar_J16_raw.index,aircraft)
fb.get_thrust_curve(polar_J16_uncorr,aircraft)
fb.correct_thrust(polar_J16_uncorr, polar_tailOn_uncorr, polar_tailOff_uncorr, dCD_CL2_tail, aircraft)
J16alpha = dir+'CT_alpha_J16.csv'
polar_J16_CT_alpha = copy.deepcopy(polar_J16_uncorr)
fb.correct_CT_alpha(polar_J16_CT_alpha,J16alpha,aircraft)
polar_J16_uncorr = copy.deepcopy(polar_J16_CT_alpha)
fb.get_uncorrect_powered(polar_J16_uncorr,aircraft)

polar_J16_corr = copy.deepcopy(polar_J16_uncorr)
Tc = np.array([point.Tc for point in polar_J16_corr.points])
fb.correct_powered(polar_J16_corr,polar_tailOff_uncorr,aircraft,tunnel)


CLu_J16 = fb.get_listfrompolar(polar_J16_uncorr,'CL')
CDu_J16 = fb.get_listfrompolar(polar_J16_uncorr,'CD')
CMu_J16 = fb.get_listfrompolar(polar_J16_uncorr,'CM')
CLc_J16 = fb.get_listfrompolar(polar_J16_corr,'CL')
CDc_J16 = fb.get_listfrompolar(polar_J16_corr,'CD')
CMc_J16 = fb.get_listfrompolar(polar_J16_corr,'CM')
CLtc_J16 = fb.get_listfrompolar(polar_J16_corr,'CLt')
CDtc_J16 = fb.get_listfrompolar(polar_J16_corr,'CDt')

CLtc_J16[:,1] = CLtc_J16[:,1]*(aircraft.wing.S/aircraft.tailh.S)

'J2.0 Case'
polar_J20 = dir+'raw_J20_PropOn_bal.txt'
polar_J20_0 = dir+'zero_PropOn_bal.txt'
polar_J20_raw = source.BAL_class.BAL_data(polar_J20,polar_J20_0,5)
polar_J20_uncorr = cP.Polar(polar_J20_raw.data,polar_J20_raw.index,aircraft)
fb.get_thrust_curve(polar_J20_uncorr,aircraft)
fb.correct_thrust(polar_J20_uncorr, polar_tailOn_uncorr, polar_tailOff_uncorr, dCD_CL2_tail, aircraft)
J20alpha = dir+'CT_alpha_J20.csv'
polar_J20_CT_alpha = copy.deepcopy(polar_J20_uncorr)
fb.correct_CT_alpha(polar_J20_CT_alpha,J20alpha,aircraft)
polar_J20_uncorr = copy.deepcopy(polar_J20_CT_alpha)
fb.get_uncorrect_powered(polar_J20_uncorr,aircraft)

polar_J20_corr = copy.deepcopy(polar_J20_uncorr)
Tc = np.array([point.Tc for point in polar_J20_corr.points])
fb.correct_powered(polar_J20_corr,polar_tailOff_uncorr,aircraft,tunnel)


CLu_J20 = fb.get_listfrompolar(polar_J20_uncorr,'CL')
CDu_J20 = fb.get_listfrompolar(polar_J20_uncorr,'CD')
CMu_J20 = fb.get_listfrompolar(polar_J20_uncorr,'CM')
CLc_J20 = fb.get_listfrompolar(polar_J20_corr,'CL')
CDc_J20 = fb.get_listfrompolar(polar_J20_corr,'CD')
CMc_J20 = fb.get_listfrompolar(polar_J20_corr,'CM')

CLtc_J20 = fb.get_listfrompolar(polar_J20_corr,'CLt')
CDtc_J20 = fb.get_listfrompolar(polar_J20_corr,'CDt')

CLtc_J20[:,1] = CLtc_J20[:,1]*(aircraft.wing.S/aircraft.tailh.S)

'J2.4 Case'
polar_J24 = dir+'raw_J24_PropOn_bal.txt'
polar_J24_0 = dir+'zero_PropOn_bal.txt'
polar_J24_raw = source.BAL_class.BAL_data(polar_J24,polar_J24_0,5)
polar_J24_uncorr = cP.Polar(polar_J24_raw.data,polar_J24_raw.index,aircraft)
fb.get_thrust_curve(polar_J24_uncorr,aircraft)
fb.correct_thrust(polar_J24_uncorr, polar_tailOn_uncorr, polar_tailOff_uncorr, dCD_CL2_tail, aircraft)
J24alpha = dir+'CT_alpha_J24.csv'
polar_J24_CT_alpha = copy.deepcopy(polar_J24_uncorr)
fb.correct_CT_alpha(polar_J24_CT_alpha,J24alpha,aircraft)
polar_J24_uncorr = copy.deepcopy(polar_J24_CT_alpha)
fb.get_uncorrect_powered(polar_J24_uncorr,aircraft)

polar_J24_corr = copy.deepcopy(polar_J24_uncorr)
Tc = np.array([point.Tc for point in polar_J24_corr.points])
fb.correct_powered(polar_J24_corr,polar_tailOff_uncorr,aircraft,tunnel)


CLu_J24 = fb.get_listfrompolar(polar_J24_uncorr,'CL')
CDu_J24 = fb.get_listfrompolar(polar_J24_uncorr,'CD')
CMu_J24 = fb.get_listfrompolar(polar_J24_uncorr,'CM')
CLc_J24 = fb.get_listfrompolar(polar_J24_corr,'CL')
CDc_J24 = fb.get_listfrompolar(polar_J24_corr,'CD')
CMc_J24 = fb.get_listfrompolar(polar_J24_corr,'CM')

CLtc_J24 = fb.get_listfrompolar(polar_J24_corr,'CLt')
CDtc_J24 = fb.get_listfrompolar(polar_J24_corr,'CDt')

CLtc_J24[:,1] = CLtc_J24[:,1]*(aircraft.wing.S/aircraft.tailh.S)

'Deflection case'
polar_delta0 = dir+'raw_delta0_bal.txt'
polar_delta0_0 = dir+'zero_tailOn_bal.txt'
polar_delta0_raw = source.BAL_class.BAL_data(polar_delta0,polar_delta0_0,5)
polar_delta0_uncorr = cP.Polar(polar_delta0_raw.data,polar_delta0_raw.index,aircraft)
fb.interpolateTc(polar_delta0_uncorr)
fb.get_uncorrect_powered(polar_delta0_uncorr, aircraft)
polar_delta0_corr = copy.deepcopy(polar_delta0_uncorr)
fb.correct_deflection(polar_delta0_corr,polar_tailOff_uncorr,polar_J16_uncorr,
                      polar_J20_uncorr,polar_J24_uncorr,polar_J16_corr,
                      polar_J20_corr,polar_J24_corr,aircraft,tunnel)


'Deflection case'
polar_delta5neg = dir+'raw_delta5neg_bal.txt'
polar_delta5neg_0 = dir+'zero_tailOn_bal.txt'
polar_delta5neg_raw = source.BAL_class.BAL_data(polar_delta5neg,polar_delta5neg_0,5)
polar_delta5neg_uncorr = cP.Polar(polar_delta5neg_raw.data,polar_delta5neg_raw.index,aircraft)
fb.interpolateTc(polar_delta5neg_uncorr)
fb.get_uncorrect_powered(polar_delta5neg_uncorr, aircraft)
polar_delta5neg_corr = copy.deepcopy(polar_delta5neg_uncorr)
fb.correct_deflection(polar_delta5neg_corr,polar_tailOff_uncorr,polar_J16_uncorr,
                      polar_J20_uncorr,polar_J24_uncorr,polar_J16_corr,
                      polar_J20_corr,polar_J24_corr,aircraft,tunnel)

'Deflection case'
polar_delta10neg = dir+'raw_delta10neg_bal.txt'
polar_delta10neg_0 = dir+'zero_tailOn_bal.txt'
polar_delta10neg_raw = source.BAL_class.BAL_data(polar_delta10neg,polar_delta10neg_0,5)
polar_delta10neg_uncorr = cP.Polar(polar_delta10neg_raw.data,polar_delta10neg_raw.index,aircraft)
fb.interpolateTc(polar_delta10neg_uncorr)
fb.get_uncorrect_powered(polar_delta10neg_uncorr, aircraft)
polar_delta10neg_corr = copy.deepcopy(polar_delta10neg_uncorr)
fb.correct_deflection(polar_delta10neg_corr,polar_tailOff_uncorr,polar_J16_uncorr,
                      polar_J20_uncorr,polar_J24_uncorr,polar_J16_corr,
                      polar_J20_corr,polar_J24_corr,aircraft,tunnel)



'Deflection case'
polar_delta5pos = dir+'raw_delta5pos_bal.txt'
polar_delta5pos_0 = dir+'zero_tailOn_bal.txt'
polar_delta5pos_raw = source.BAL_class.BAL_data(polar_delta5pos,polar_delta5pos_0,5)
polar_delta5pos_uncorr = cP.Polar(polar_delta5pos_raw.data,polar_delta5pos_raw.index,aircraft)
fb.interpolateTc(polar_delta5pos_uncorr)
fb.get_uncorrect_powered(polar_delta5pos_uncorr, aircraft)
polar_delta5pos_corr = copy.deepcopy(polar_delta5pos_uncorr)
fb.correct_deflection(polar_delta5pos_corr,polar_tailOff_uncorr,polar_J16_uncorr,
                      polar_J20_uncorr,polar_J24_uncorr,polar_J16_corr,
                      polar_J20_corr,polar_J24_corr,aircraft,tunnel)
'Deflection case'
polar_delta10pos = dir+'raw_delta10pos_bal.txt'
polar_delta10pos_0 = dir+'zero_tailOn_bal.txt'
polar_delta10pos_raw = source.BAL_class.BAL_data(polar_delta10pos,polar_delta10pos_0,5)
polar_delta10pos_uncorr = cP.Polar(polar_delta10pos_raw.data,polar_delta10pos_raw.index,aircraft)
fb.interpolateTc(polar_delta10pos_uncorr)
fb.get_uncorrect_powered(polar_delta10pos_uncorr, aircraft)
polar_delta10pos_corr = copy.deepcopy(polar_delta10pos_uncorr)
fb.correct_deflection(polar_delta10pos_corr,polar_tailOff_uncorr,polar_J16_uncorr,
                      polar_J20_uncorr,polar_J24_uncorr,polar_J16_corr,
                      polar_J20_corr,polar_J24_corr,aircraft,tunnel)

CMc_delta0 = np.array([point.CMp for point in polar_delta0_corr.points])
CMc_delta10neg = np.array([point.CMp for point in polar_delta10neg_corr.points])
CMc_delta5neg = np.array([point.CMp for point in polar_delta5neg_corr.points])
CMc_delta5pos = np.array([point.CMp for point in polar_delta5pos_corr.points])
CMc_delta10pos = np.array([point.CMp for point in polar_delta10pos_corr.points])

CLc_delta0 = np.array([point.CFl for point in polar_delta0_corr.points])
CLc_delta10neg = np.array([point.CFl for point in polar_delta10neg_corr.points])
CLc_delta5neg = np.array([point.CFl for point in polar_delta5neg_corr.points])
CLc_delta5pos = np.array([point.CFl for point in polar_delta5pos_corr.points])
CLc_delta10pos = np.array([point.CFl for point in polar_delta10pos_corr.points])

CMc_cg_delta0 = CMc_delta0 + CLc_delta0*((0.5-0.25)*aircraft.wing.cmac)
CMc_cg_delta10neg = CMc_delta10neg + CLc_delta10neg*((0.5-0.25)*aircraft.wing.cmac)
CMc_cg_delta5neg = CMc_delta5neg + CLc_delta5neg*((0.5-0.25)*aircraft.wing.cmac)
CMc_cg_delta5pos = CMc_delta5pos + CLc_delta5pos*((0.5-0.25)*aircraft.wing.cmac)
CMc_cg_delta10pos = CMc_delta10pos + CLc_delta10pos*((0.5-0.25)*aircraft.wing.cmac)

delta_e = np.array([-10,-5,0,5,10])
CLc_delta_J16 = np.array([CLc_delta10neg[0],CLc_delta5neg[0],CLc_delta0[0],
                          CLc_delta5pos[0],CLc_delta10pos[0]])
CLc_delta_J20 = np.array([CLc_delta10neg[1],CLc_delta5neg[1],CLc_delta0[1],
                          CLc_delta5pos[1],CLc_delta10pos[1]])
CLc_delta_J24 = np.array([CLc_delta10neg[2],CLc_delta5neg[2],CLc_delta0[2],
                          CLc_delta5pos[2],CLc_delta10pos[2]])

CMc_cg_delta_J16 = np.array([CMc_cg_delta10neg[0],CMc_cg_delta5neg[0],CMc_cg_delta0[0],
                          CMc_cg_delta5pos[0],CMc_cg_delta10pos[0]])
CMc_cg_delta_J20 = np.array([CMc_cg_delta10neg[1],CMc_cg_delta5neg[1],CMc_cg_delta0[1],
                          CMc_cg_delta5pos[1],CMc_cg_delta10pos[1]])
CMc_cg_delta_J24 = np.array([CMc_cg_delta10neg[2],CMc_cg_delta5neg[2],CMc_cg_delta0[2],
                          CMc_cg_delta5pos[2],CMc_cg_delta10pos[2]])




plt.rcParams.update({'font.size':14})


plt.figure(1)
plt.plot(delta_e,CLc_delta_J16,'--x')
plt.plot(delta_e,CLc_delta_J20,'--x')
plt.plot(delta_e,CLc_delta_J24,'--x')
plt.ylabel(r'Lift coefficient $C_{L}$')
plt.xlabel(r'Elevator deflection angle $\delta_{e}$ [deg]')
plt.legend([r'$C_{L,J=1.6}$',r'$C_{L,J=2.0}$',r'$C_{L,J=2.4}$'])
plt.grid(b=True, which='major',linestyle='-')
plt.grid(b=True, which='minor',linestyle='--')
plt.minorticks_on()

plt.figure(2)
plt.plot(delta_e,CMc_cg_delta_J16,'--x')
plt.plot(delta_e,CMc_cg_delta_J20,'--x')
plt.plot(delta_e,CMc_cg_delta_J24,'--x')
plt.legend([r'$C_{M,J=1.6}$',r'$C_{M,J=2.0}$',r'$C_{M,J=2.4}$'])
plt.ylabel(r'Moment coefficient about c.g. $C_{M,c.g.}$')
plt.xlabel(r'Elevator deflection angle $\delta_{e}$ [deg]')
plt.legend([r'$C_{M,J=1.6}$',r'$C_{M,J=2.0}$',r'$C_{M,J=2.4}$'])
plt.grid(b=True, which='major',linestyle='-')
plt.grid(b=True, which='minor',linestyle='--')
plt.minorticks_on()