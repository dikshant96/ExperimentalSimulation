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
CLtc_J16_0 = np.array([point.CFl_t_0 for point in polar_J16_corr.points])

CLtc_J16[:,1] = CLtc_J16[:,1]*(aircraft.wing.S/aircraft.tailh.S)
CLtc_J16_0 = CLtc_J16_0*(aircraft.wing.S/aircraft.tailh.S)

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
CLtc_J20_0 = np.array([point.CFl_t_0 for point in polar_J20_corr.points])


CLtc_J20[:,1] = CLtc_J20[:,1]*(aircraft.wing.S/aircraft.tailh.S)
CLtc_J20_0 = CLtc_J20_0*(aircraft.wing.S/aircraft.tailh.S)


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
CLtc_J24_0 = np.array([point.CFl_t_0 for point in polar_J24_corr.points])

CLtc_J24[:,1] = CLtc_J24[:,1]*(aircraft.wing.S/aircraft.tailh.S)
CLtc_J24_0 = CLtc_J24_0*(aircraft.wing.S/aircraft.tailh.S)


'Translate moments to cg'
CMc_cg_tailOn = copy.deepcopy(CMc_tailOn)
CMc_cg_tailOn[:,1] = CMc_tailOn[:,1] + CLc_tailOn[:,1]*((0.5-0.25)*aircraft.wing.cmac)

CMc_cg_J16 = copy.deepcopy(CMc_J16)
Fz_J16 = np.array([point.F3 for point in polar_J16_corr.points])
qInf_J16 = np.array([point.qInf for point in polar_J16_corr.points])
CFz_J16 = Fz_J16/(qInf_J16*aircraft.wing.S)
CMc_cg_J16[:,1] = CMc_J16[:,1] + CFz_J16*((0.5-0.25)*aircraft.wing.cmac)

CMc_cg_J20 = copy.deepcopy(CMc_J20)
Fz_J20 = np.array([point.F3 for point in polar_J20_corr.points])
qInf_J20 = np.array([point.qInf for point in polar_J20_corr.points])
CFz_J20 = Fz_J20/(qInf_J20*aircraft.wing.S)
CMc_cg_J20[:,1] = CMc_J20[:,1] + CFz_J20*((0.5-0.25)*aircraft.wing.cmac)

CMc_cg_J24 = copy.deepcopy(CMc_J24)
Fz_J24 = np.array([point.F3 for point in polar_J24_corr.points])
qInf_J24 = np.array([point.qInf for point in polar_J24_corr.points])
CFz_J24 = Fz_J24/(qInf_J24*aircraft.wing.S)
CMc_cg_J24[:,1] = CMc_J24[:,1] + CFz_J24*((0.5-0.25)*aircraft.wing.cmac)

plt.rcParams.update({'font.size':14})

# plt.figure(1)
# plt.plot(CMc_cg_tailOn[:-3,0],CMc_cg_tailOn[:-3,1],'--x')
# plt.plot(CMc_cg_J16[:,0],CMc_cg_J16[:,1],'--x')
# plt.plot(CMc_cg_J20[:,0],CMc_cg_J20[:,1],'--x')
# plt.plot(CMc_cg_J24[:,0],CMc_cg_J24[:,1],'--x')
# plt.legend([r'$C_{M,propOff}$',r'$C_{M,J=1.6}$',r'$C_{M,J=2.0}$',r'$C_{M,J=2.4}$'])
# plt.ylabel(r'Moment coefficient about c.g. $C_{M,c.g.}$')
# plt.xlabel(r'Angle of attack $\alpha$ [deg]')
# plt.grid(b=True, which='major',linestyle='-')
# plt.grid(b=True, which='minor',linestyle='--')
# plt.minorticks_on()
#
# plt.figure(2)
# plt.plot(CLtc_tailOn[:-3,0],CLtc_tailOn[:-3,1],'--x')
# plt.plot(CLtc_J16[:,0],CLtc_J16[:,1],'--x')
# plt.plot(CLtc_J20[:,0],CLtc_J20[:,1],'--x')
# plt.plot(CLtc_J24[:,0],CLtc_J24[:,1],'--x')
# plt.legend([r'$C_{L,t,propOff}$',r'$C_{L,t,J=1.6}$',r'$C_{L,t,J=2.0}$',r'$C_{L,t,J=2.4}$'])
# plt.ylabel(r'Lift coefficient $C_{L}$')
# plt.xlabel(r'Angle of attack $\alpha$ [deg]')
# plt.grid(b=True, which='major',linestyle='-')
# plt.grid(b=True, which='minor',linestyle='--')
# plt.minorticks_on()
#
# plt.figure(3)
# plt.plot(CLtc_tailOn[:-3,0],CLtc_tailOn[:-3,1],'--x')
# plt.plot(CLtc_J16[:,0],CLtc_J16_0,'--x')
# plt.plot(CLtc_J20[:,0],CLtc_J20_0,'--x')
# plt.plot(CLtc_J24[:,0],CLtc_J24_0,'--x')
# plt.legend([r'$C_{L,t,propOff}$',r'$C_{L,t,J=1.6}$',r'$C_{L,t,J=2.0}$',r'$C_{L,t,J=2.4}$'])
# plt.ylabel(r'Lift coefficient $C_{L}$')
# plt.xlabel(r'Angle of attack $\alpha$ [deg]')
# plt.grid(b=True, which='major',linestyle='-')
# plt.grid(b=True, which='minor',linestyle='--')
# plt.minorticks_on()
# plt.figure(1)
# plt.plot(CMc_J16[:,0],CMc_J16[:,1],'--x')
# plt.plot(CMu_J16[:,0],CMu_J16[:,1],'--x')
# plt.legend([r'$C_{M,c}$',r'$C_{M,u}$'])
# plt.ylabel(r'Moment coefficient about c/4 $C_{M,c/4}$')
# plt.xlabel(r'Angle of attack $\alpha$ [deg]')
# plt.grid(b=True, which='major',linestyle='-')
# plt.grid(b=True, which='minor',linestyle='--')
# plt.minorticks_on()
#
# plt.figure(2)
# plt.plot(CLc_J16[:,0],CLc_J16[:,1],'--x')
# plt.plot(CLu_J16[:,0],CLu_J16[:,1],'--o')
# plt.legend([r'$C_{L,c}$',r'$C_{L,u}$'])
# plt.ylabel(r'Lift coefficient $C_{L}$')
# plt.xlabel(r'Angle of attack $\alpha$ [deg]')
# plt.grid(b=True, which='major',linestyle='-')
# plt.grid(b=True, which='minor',linestyle='--')
# plt.minorticks_on()
#
# plt.figure(3)
# plt.plot(CDc_J16[:,0],CDc_J16[:,1],'--x')
# plt.plot(CDu_J16[:,0],CDu_J16[:,1],'--o')
# #plt.plot(CDtc_tailOn[:,0],CDtc_tailOn[:,1]*(aircraft.wing.S/aircraft.tailh.S),'--s')
# plt.legend([r'$C_{D,c}$',r'$C_{D,u}$'])
# plt.ylabel(r'Drag coefficient $C_{D}$')
# plt.xlabel(r'Angle of attack $\alpha$ [deg]')
# plt.grid(b=True, which='major',linestyle='-')
# plt.grid(b=True, which='minor',linestyle='--')
# plt.minorticks_on()