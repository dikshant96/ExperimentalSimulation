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

'Shift moment to c.g.'

CMc_cg_tailOn = copy.deepcopy(CMc_tailOn)
CMc_cg_tailOn[:,1] = CMc_tailOn[:,1] + CLc_tailOn[:,1]*((0.5-0.25)*aircraft.wing.cmac)
dCMc_cg_alpha, c = np.polyfit(CMc_cg_tailOn[:,0],CMc_cg_tailOn[:,1],1)
dCMc_cg_alpha_2, c2 = np.polyfit(CMc_cg_tailOn[:-3,0],CMc_cg_tailOn[:-3,1],1)
CMc_cg_tailOn_fit = dCMc_cg_alpha*CMc_cg_tailOn[:,0] + c
CMc_cg_tailOn_fit_2 = dCMc_cg_alpha_2*CMc_cg_tailOn[:,0] + c2

plt.rcParams.update({'font.size':14})

# plt.figure(1)
# plt.plot(CMc_cg_tailOn[:,0],CMc_cg_tailOn[:,1],'--x')
# plt.plot(CMc_cg_tailOn[:,0],CMc_cg_tailOn_fit[:],'-')
# plt.plot(CMc_cg_tailOn[:,0],CMc_cg_tailOn_fit_2[:],'-')
# plt.legend([r'$C_{M}$',r'$C_{M,fit}$',r'$C_{M,nostall}$'])
# plt.ylabel(r'Moment coefficient about c.g. $C_{M,c.g.}$')
# plt.xlabel(r'Angle of attack $\alpha$ [deg]')
# plt.grid(b=True, which='major',linestyle='-')
# plt.grid(b=True, which='minor',linestyle='--')
# plt.minorticks_on()
#
# plt.figure(2)
# plt.plot(CLc_tailOn[:,0],CLc_tailOn[:,1],'--x')
# plt.plot(CLwc_tailOn[:,0],CLwc_tailOn[:,1],'--o')
# plt.plot(CLtc_tailOn[:,0],CLtc_tailOn[:,1]*(aircraft.wing.S/aircraft.tailh.S),'--s')
# plt.legend([r'$C_{L}$',r'$C_{L,w}$',r'$C_{L,t}$'])
# plt.ylabel(r'Lift coefficient $C_{L}$')
# plt.xlabel(r'Angle of attack $\alpha$ [deg]')
# plt.grid(b=True, which='major',linestyle='-')
# plt.grid(b=True, which='minor',linestyle='--')
# plt.minorticks_on()
#
# plt.figure(3)
# plt.plot(CDc_tailOn[:,0],CDc_tailOn[:,1],'--x')
# plt.plot(CDwc_tailOn[:,0],CDwc_tailOn[:,1],'--o')
# plt.plot(CDtc_tailOn[:,0],CDtc_tailOn[:,1]*(aircraft.wing.S/aircraft.tailh.S),'--s')
# plt.legend([r'$C_{D}$',r'$C_{D,w}$',r'$C_{D,t}$'])
# plt.ylabel(r'Drag coefficient $C_{D}$')
# plt.xlabel(r'Angle of attack $\alpha$ [deg]')
# plt.grid(b=True, which='major',linestyle='-')
# plt.grid(b=True, which='minor',linestyle='--')
# plt.minorticks_on()

plt.figure(1)
plt.plot(CMc_tailOff[:,0],CMc_tailOff[:,1],'--x')
plt.plot(CMu_tailOff[:,0],CMu_tailOff[:,1],'--x')
plt.legend([r'$C_{M,c}$',r'$C_{M,u}$'])
plt.ylabel(r'Moment coefficient about c/4 $C_{M,c/4}$')
plt.xlabel(r'Angle of attack $\alpha$ [deg]')
plt.grid(b=True, which='major',linestyle='-')
plt.grid(b=True, which='minor',linestyle='--')
plt.minorticks_on()

plt.figure(2)
plt.plot(CLc_tailOff[:,0],CLc_tailOff[:,1],'--x')
plt.plot(CLu_tailOff[:,0],CLu_tailOff[:,1],'--o')
plt.legend([r'$C_{L,c}$',r'$C_{L,u}$'])
plt.ylabel(r'Lift coefficient $C_{L}$')
plt.xlabel(r'Angle of attack $\alpha$ [deg]')
plt.grid(b=True, which='major',linestyle='-')
plt.grid(b=True, which='minor',linestyle='--')
plt.minorticks_on()

plt.figure(3)
plt.plot(CDc_tailOff[:,0],CDc_tailOff[:,1],'--x')
plt.plot(CDu_tailOff[:,0],CDu_tailOff[:,1],'--o')
#plt.plot(CDtc_tailOff[:,0],CDtc_tailOff[:,1]*(aircraft.wing.S/aircraft.tailh.S),'--s')
plt.legend([r'$C_{D,c}$',r'$C_{D,u}$'])
plt.ylabel(r'Drag coefficient $C_{D}$')
plt.xlabel(r'Angle of attack $\alpha$ [deg]')
plt.grid(b=True, which='major',linestyle='-')
plt.grid(b=True, which='minor',linestyle='--')
plt.minorticks_on()