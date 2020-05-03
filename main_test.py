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
polar_tailOff_corr_block = cP.Polar(raw_data_40,polar_tailOff_raw.index,aircraft)
fb.correct_blockage(polar_tailOff_corr_block,aircraft,tunnel,'w')
polar_tailOff_corr = copy.deepcopy(polar_tailOff_corr_block)
fb.correct_streamline(polar_tailOff_corr,aircraft,tunnel,'w')

df = pd.DataFrame(np.zeros((len(polar_tailOff_corr.points),11)),
                  columns=['alphau','alpha','CLu','CDu','CMu',
                          'CLc_b','CDc_b','CMc_b','CL','CD','CM'])
for i, point in enumerate(polar_tailOff_corr.points):
    df.at[i,'alphau'] = polar_tailOff_uncorr.points[i].alpha
    df.at[i,'alpha'] = polar_tailOff_corr.points[i].alpha
    df.at[i,'CLu'] = polar_tailOff_uncorr.points[i].CFl
    df.at[i,'CDu'] = polar_tailOff_uncorr.points[i].CFd
    df.at[i,'CMu'] = polar_tailOff_uncorr.points[i].CMp
    df.at[i,'CLc_b'] = polar_tailOff_corr_block.points[i].CFl
    df.at[i,'CDc_b'] = polar_tailOff_corr_block.points[i].CFd
    df.at[i,'CMc_b'] = polar_tailOff_corr_block.points[i].CMp
    df.at[i,'CL'] = polar_tailOff_corr.points[i].CFl
    df.at[i,'CD'] = polar_tailOff_corr.points[i].CFd
    df.at[i,'CM'] = polar_tailOff_corr.points[i].CMp

plt.figure(1)
plt.plot(df[:]['alphau'],df[:]['CLu'],'--')
plt.plot(df[:]['alphau'],df[:]['CLc_b'],'--')
plt.plot(df[:]['alpha'],df[:]['CL'],'--')

plt.figure(2)
plt.plot(df[:]['alphau'],df[:]['CDu'],'--')
plt.plot(df[:]['alphau'],df[:]['CDc_b'],'--')
plt.plot(df[:]['alpha'],df[:]['CD'],'--')

plt.figure(3)
plt.plot(df[:]['alphau'],df[:]['CMu'],'--')
plt.plot(df[:]['alphau'],df[:]['CMc_b'],'--')
plt.plot(df[:]['alpha'],df[:]['CM'],'--')