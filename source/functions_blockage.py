import numpy as np
import matplotlib.pyplot as plt

def get_solidblockage(aircraft,tunnel,bool):
    t1w = 0.87
    k1w = 1.02
    esb_w = (k1w * t1w * aircraft.wing.V) / ((tunnel.C) ** (3 / 2.))

    t1f = 0.86
    k3f = 0.915
    esb_f = (k3f * t1f * aircraft.fuselage.V) / ((tunnel.C) ** (3 / 2.))

    t1ss = 0.86
    k1ss = 0.90
    ss_V = 0.0035296
    esb_ss = (k1ss * t1ss * ss_V) / ((tunnel.C) ** (3 / 2.))

    t1ms = 0.86
    k1ms = 0.90
    ms_V = 0.0004491
    esb_ms = (t1ms * k1ms * ms_V) / ((tunnel.C) ** (3 / 2.))

    t1h = 0.86
    k1h = 1.035
    esb_th = (t1h * k1h * aircraft.tailh.V) / ((tunnel.C) ** (3 / 2.))

    t1v = 0.86
    k1v = 1.035
    esb_tv = (t1v * k1v * aircraft.tailv.V) / ((tunnel.C) ** (3 / 2.))

    if bool == 'w':
        esb = esb_f + esb_w + esb_ms + esb_ss
    elif bool == 'wt':
        esb = esb_f + esb_w + esb_ss + esb_ms + esb_th + esb_tv

    return esb

def get_wakeblockage(polar,aircraft,tunnel):
    CLu = np.array([point.CFl for point in polar.points])
    CDu = np.array([point.CFd for point in polar.points])
    CD0 = np.min(CDu)
    CLu_2 = CLu**2
    CLu_2_linear = []
    CDu_linear = []
    for i,CL in enumerate(CLu_2):
        if 0.0 < CL < 0.7:
            CLu_2_linear.append(CL)
            CDu_linear.append(CDu[i])
    m,c = np.polyfit(CLu_2_linear,CDu_linear,1)
    CDi = m*CLu_2
    CDs_t = CDu - CD0 - CDi
    CDs = np.array([max(CD, 0) for CD in CDs_t])

    ewb = (aircraft.wing.S/(4*tunnel.C)) * CD0 + ((5 * aircraft.wing.S)/(4 * tunnel.C)) * CDs
    return ewb

def correct_blockage(polar,aircraft,tunnel,bool):
    esb = get_solidblockage(aircraft, tunnel, bool)
    ewb = get_wakeblockage(polar,aircraft, tunnel)
    for i,point in enumerate(polar.points):
        point.et = esb + ewb[i]
        point.qInf = point.qInf*((1+point.et)**2)
        point.U = point.U*(1+point.et)
        point.get_coeffs(aircraft)
        point.get_modelaxiscoeffs(aircraft)
        point.get_conventionalcoeffs(aircraft)

def get_streamline_wing(CLu,alphau,aircraft,tunnel):
    bv = 0.76 * aircraft.wing.b
    be = (aircraft.wing.b + bv) / 2
    delta = 0.108
    delta_alpha_wing_sc = delta * (aircraft.wing.S / tunnel.C) * (180 / np.pi) * CLu
    t2w = 0.11
    delta_alpha_wing = (1 + t2w) * delta_alpha_wing_sc
    delta_drag_wing = delta * (aircraft.wing.S / tunnel.C) * CLu ** 2
    CLu_linear = []
    alphau_linear = []
    for i, alpha_i in enumerate(alphau):
        if -1 < alpha_i < 7:
            alphau_linear.append(alpha_i)
            CLu_linear.append(CLu[i])
    CLu_linear = np.array(CLu_linear)
    alphau_linear = np.array(alphau_linear)
    CLalpha, c = np.polyfit(alphau_linear, CLu_linear, 1)
    delta_moment_wing = 0.125 * delta_alpha_wing * CLalpha
    delta_array_wing = np.array([[delta_alpha_wing[i],delta_drag_wing[i],delta_moment_wing[i]] for i in range(len(delta_alpha_wing))])
    return delta_array_wing

def correct_notail(polar,aircraft,tunnel,bool):
    esb = get_solidblockage(aircraft,tunnel, bool)
    ewb = get_wakeblockage(polar,aircraft, tunnel)
    CLu = np.array([point.CFl for point in polar.points])
    alphau = np.array([point.alpha for point in polar.points])
    delta_array_wing = get_streamline_wing(CLu, alphau, aircraft, tunnel)
    for i,point in enumerate(polar.points):
        point.et = esb + ewb[i]
        point.qInf = point.qInf * ((1 + point.et) ** 2)
        point.U = point.U * (1 + point.et)
        point.get_reynolds(aircraft)
        point.alpha = point.alpha + delta_array_wing[i, 0]
        point.CFd = point.CFd*(1/(1+point.et)**2) + delta_array_wing[i, 1]
        point.CMp = point.CMp*(1/(1+point.et)**2) + delta_array_wing[i, 2]
        point.CFl = point.CFl*(1/(1+point.et)**2)

def get_streamline_tail(CMu,CMwu,CLwu,alphau,aircraft,tunnel):
    CMh = CMu - CMwu
    CMh_linear = []
    alpha_linear = []
    for i, alpha in enumerate(alphau):
        if -1 < alpha < 7:
            alpha_linear.append(alpha)
            CMh_linear.append(CMh[i])
    CMhalpha, c = np.polyfit(alpha_linear,CMh_linear,1)
    #print(CMhalpha)
    lth = 3.22*aircraft.wing.cmac
    delta = 0.108
    t2h = 0.88
    #print(CLwu)
    delta_alpha_tail = delta*t2h*(aircraft.wing.S/tunnel.C)*CLwu*(180/np.pi)
    delta_moment_tail = delta_alpha_tail*CMhalpha
    delta_array_tail = np.array([[delta_alpha_tail[i], delta_moment_tail[i]] for i in range(len(delta_alpha_tail))])
    return delta_array_tail

def correct_tail(polar,polar_tailOff_uncorr,aircraft,tunnel,bool):
    esb = get_solidblockage(aircraft,tunnel,'wt')
    ewb = get_wakeblockage(polar,aircraft,tunnel)

    'Uncorrected tail off'
    CLwu = np.array([point.CFl for point in polar_tailOff_uncorr.points])
    alphawu = np.array([point.alpha for point in polar_tailOff_uncorr.points])
    CMwu = np.array([point.CMp for point in polar_tailOff_uncorr.points])
    CDwu = np.array([point.CFd for point in polar_tailOff_uncorr.points])
    delta_array_wing = get_streamline_wing(CLwu, alphawu, aircraft, tunnel)

    'Uncorrected net'
    CLu = np.array([point.CFl for point in polar.points])
    CDu = np.array([point.CFd for point in polar.points])
    CMu = np.array([point.CMp for point in polar.points])
    alphau = np.array([point.alpha for point in polar.points])

    'Get same alphas for uncorrected wing'
    CMwu_t = []
    CLwu_t = []
    alphawu_t = []
    CDwu_t = []
    delta_array_wing_tail = []
    for i, alpha in enumerate(alphau):
        min_index = np.argmin(np.abs(alphawu - alpha))
        CMwu_t.append(CMwu[min_index])
        CLwu_t.append(CLwu[min_index])
        alphawu_t.append(alphawu[min_index])
        CDwu_t.append(CDwu[min_index])
        delta_array_wing_tail.append(delta_array_wing[min_index,:])

    CMwu_t = np.array(CMwu_t)
    CLwu_t = np.array(CLwu_t)
    CDwu_t = np.array(CDwu_t)
    alphawu_t = np.array(alphawu_t)
    delta_array_wing_tail = np.array(delta_array_wing_tail)

    delta_array_tail  = get_streamline_tail(CMu,CMwu_t,CLwu_t,alphau,aircraft,tunnel)
    for i,point in enumerate(polar.points):
        point.et = esb + ewb[i]
        point.qInf = point.qInf * ((1 + point.et) ** 2)
        point.U = point.U * (1 + point.et)
        point.get_reynolds(aircraft)
        alpha_temp = point.alpha
        point.alpha = point.alpha + delta_array_wing_tail[i, 0]
        point.alpha_t = alpha_temp + delta_array_tail[i, 0]
        point.CFL_w = CLwu_t[i] * (1 / (1 + point.et) ** 2)
        point.CFd_w = CDwu_t[i] * (1 / (1 + point.et) ** 2) + delta_array_wing_tail[i, 1]
        point.CMp_w = CMwu_t[i] * (1 / (1 + point.et) ** 2) + delta_array_wing_tail[i, 2]
        point.CFd = point.CFd * (1 / (1 + point.et) ** 2) + delta_array_wing_tail[i, 1]
        point.CMp = point.CMp * (1 / (1 + point.et) ** 2) + delta_array_wing_tail[i, 2] + delta_array_tail[i, 1]
        point.CFl = point.CFl * (1 / (1 + point.et) ** 2)
        point.CFl_t = point.CFl - point.CFL_w
        point.CFd_t = point.CFd - point.CFd_w
        point.CMp_t = point.CMp - point.CMp_w

def get_listfrompolar(polar,bool):
    alpha = []
    var = []
    for i, point in enumerate(polar.points):
        if bool == 'CL':
            alpha.append(point.alpha)
            var.append(point.CFl)
        elif bool == 'CD':
            alpha.append(point.alpha)
            var.append(point.CFd)
        elif bool == 'CM':
            alpha.append(point.alpha)
            var.append(point.CMp)
        elif bool == 'CLw':
            alpha.append(point.alpha)
            var.append(point.CFL_w)
        elif bool == 'CDw':
            alpha.append(point.alpha)
            var.append(point.CFd_w)
        elif bool == 'CMw':
            alpha.append(point.alpha)
            var.append(point.CMp_w)
        elif bool == 'CLt':
            alpha.append(point.alpha_t)
            var.append(point.CFl_t)
        elif bool == 'CDt':
            alpha.append(point.alpha_t)
            var.append(point.CFd_t)
        elif bool == 'CMt':
            alpha.append(point.alpha_t)
            var.append(point.CMp_t)
    array = np.array([[alpha[i], var[i]] for i in range(len(var))])
    return array

def get_thrust_curve(polar,aircraft):
    J = np.array([point.J for point in polar.points])
    rho = np.array([point.rho for point in polar.points])
    n = np.array([point.n for point in polar.points])
    qInf = np.array([point.qInf for point in polar.points])
    J_fit = np.array([1.756,1.873,1.982,2.073,2.194,2.299])
    CT_fit = np.array([0.236, 0.195, 0.153, 0.119, 0.075, 0.028])
    grad, c = np.polyfit(J_fit,CT_fit,1)
    CT_calc = grad*J + c
    T_calc = CT_calc*rho*(n**2)*(aircraft.prop.D**4)
    Sp = np.pi * (aircraft.prop.D/2)**2
    Tc_calc = T_calc/(qInf*Sp)
    for i, point in enumerate(polar.points):
        point.CT = CT_calc[i]
        point.T = T_calc[i]
        point.Tc = Tc_calc[i]

# def get_thrust_zero(polar,polar_propOff,aircraft):
#     index_0 = np.argmin(np.abs(np.array([point.alpha for point in polar.points])))
#     index_0_propOff = np.argmin(np.abs(np.array([point.alpha for point in polar_propOff.points])))
#     Fx = polar.points[index_0].F1
#     Fx_propOff = polar_propOff.points[index_0_propOff].F1
#     Fz = polar.points[index_0].F3
#     Fz_propOff = polar_propOff.points[index_0_propOff].F3
#     #print(index_0,index_0_propOff)
#     L = Fz_propOff
#     N = Fz - L
#     D = Fx_propOff
#     T = Fx_propOff - Fx
#     CT = (T/2)/(polar.points[index_0].rho*polar.points[index_0].n**2*aircraft.prop.D**4)
#     return L, N, D, T, CT

def get_slipstream(Tc, aircraft):
    var1 = Tc*np.sqrt((np.sqrt(1 + Tc) + 1)/(2*np.sqrt(1+Tc)))
    var2 = aircraft.prop.Nblades*0.6*(aircraft.prop.D/aircraft.tailh.b)*1
    q_ratio = 1 + var2*var1
    return q_ratio


def get_CT_alpha(fname,alpha_i):
    f = open(fname,'r')
    lines = f.readlines()
    f.close()
    CT = []
    alpha = []
    for line in lines:
        line = line.strip()
        line = line.split(',')
        alpha.append(float(line[0]))
        CT.append(float(line[1]))
    CT_alpha_ratio = np.interp(alpha_i, alpha, CT)/np.interp(0, alpha, CT)
    return CT_alpha_ratio

def correct_CT_alpha(polar,Jfname,aircraft):
    for point in polar.points:
        alpha_i = point.alpha
        CT_ratio = get_CT_alpha(Jfname,alpha_i)
        point.CT = point.CT*CT_ratio
        T_calc = point.CT * point.rho * (point.n ** 2) * (aircraft.prop.D ** 4)
        Sp = np.pi * (aircraft.prop.D / 2) ** 2
        Tc_calc = T_calc / (point.qInf * Sp)
        point.T = T_calc
        point.Tc = Tc_calc

def get_tail_CLalpha(CLtc_tailOn,aircraft):
    dCLtc, c = np.polyfit(CLtc_tailOn[:,0],CLtc_tailOn[:,1]*(aircraft.wing.S/aircraft.tailh.S),1)
    return dCLtc

def get_tail_CDCL2(CDtc_tailOn,CLtc_tailOn,aircraft):
    CD = []
    CL2 = []
    for i in range(len(CLtc_tailOn)):
        if -0.1 < CLtc_tailOn[i,1] < 0.2:
            CL2.append(CLtc_tailOn[i,1]**2)
            CD.append(CDtc_tailOn[i,1])
    #print(CL2)
    dCD_CL2, c = np.polyfit(CL2,CD,1)
    return dCD_CL2


def correct_thrust(polar,polar_tailOn,polar_tailOff,dCD_CL2_tail,aircraft):
    diff = 1
    iter = 1

    point_0 = polar.points[1]
    point_0_propOff = polar_tailOn.points[2]
    point_0_tailOff = polar_tailOff.points[5]
    Tc = point_0.Tc
    while diff > 0.001 and iter < 20:
        Fx = point_0.F1
        Dwt = point_0_propOff.F1
        Sp = np.pi * (aircraft.prop.D / 2) ** 2
        T_calc = Tc * point_0.qInf * Sp
        q_ratio_0 = get_slipstream(Tc, aircraft)
        Lwt = point_0_propOff.F3
        Lw = point_0_tailOff.F3
        Lt = Lwt - Lw
        Fz = point_0.F3
        delta_Lt = Fz - Lwt
        # delta_Lt_slipstream = Lt * (q_ratio_0 - 1)
        # delta_CLt = (delta_Lt_slipstream/(aircraft.wing.S * point_0.qInf))*(aircraft.wing.S/aircraft.tailh.S)
        delta_CLt = (delta_Lt / (aircraft.wing.S * point_0.qInf)) * (aircraft.wing.S / aircraft.tailh.S)
        delta_CDt = (delta_CLt ** 2) * dCD_CL2_tail
        delta_Dt = delta_CDt * point_0.qInf * aircraft.tailh.S
        T = Dwt + delta_Dt - Fx
        T = T / 2
        diff = np.abs(T_calc - T)
        iter = iter + 1
        Tc_old = Tc
        Tc_new = T / (point_0.qInf * Sp)
        Tc = Tc_old * 0.75 + 0.25 * Tc_new
    for point in polar.points:
        point.Tc = Tc
        point.T = Tc*(point_0.qInf*Sp)
        point.CT = point.T/(point.rho * (point.n ** 2) * (aircraft.prop.D ** 4))


def get_slipstreamblockage(polar,aircraft,tunnel):
    Tc = np.array([point.Tc for point in polar.points])
    var1 = np.sqrt(1 + (8/np.pi)*Tc)-1
    Sp = np.pi*(aircraft.prop.D/2)**2
    ess = (-1/(4*np.pi*tunnel.b*tunnel.h))*Sp*var1
    return ess

def get_uncorrected_drag(polar,aircraft):
    for point in polar.points:
        T = 2*point.T
        #T = T*1.03
        point.D = point.F1 + T*np.cos(point.alpha*(np.pi/180))
        point.CFd = point.D/(point.qInf * aircraft.wing.S)

def get_uncorrect_powered(polar,aircraft):
    Sp = np.pi * (aircraft.prop.D / 2) ** 2
    for i,point in enumerate(polar.points):
        T = point.Tc*point.qInf*Sp
        T2 = 2*T
        D = point.F1 + T2*np.cos(point.alpha*(np.pi/180))
        L = point.F3 - T2*np.sin(point.alpha*(np.pi/180))
        point.CFd = D/(point.qInf * aircraft.wing.S)
        point.CFl = L/(point.qInf * aircraft.wing.S)


def correct_powered(polar,polar_tailOff,aircraft,tunnel):
    esb = get_solidblockage(aircraft, tunnel, 'wt')
    ewb = get_wakeblockage(polar, aircraft, tunnel)
    #Tc =  np.array([point.Tc for point in polar.points])
    #print('Tc')
    ess = get_slipstreamblockage(polar, aircraft, tunnel)

    'Uncorrected tail off'
    CLwu = np.array([point.CFl for point in polar_tailOff.points])
    alphawu = np.array([point.alpha for point in polar_tailOff.points])
    CMwu = np.array([point.CMp for point in polar_tailOff.points])
    CDwu = np.array([point.CFd for point in polar_tailOff.points])
    delta_array_wing = get_streamline_wing(CLwu, alphawu, aircraft, tunnel)

    'Uncorrected net'
    CLu = np.array([point.CFl for point in polar.points])
    CDu = np.array([point.CFd for point in polar.points])
    CMu = np.array([point.CMp for point in polar.points])
    alphau = np.array([point.alpha for point in polar.points])

    'Get same alphas for uncorrected wing'
    CMwu_t = []
    CLwu_t = []
    alphawu_t = []
    CDwu_t = []
    delta_array_wing_tail = []
    for i, alpha in enumerate(alphau):
        min_index = np.argmin(np.abs(alphawu - alpha))
        CMwu_t.append(CMwu[min_index])
        CLwu_t.append(CLwu[min_index])
        alphawu_t.append(alphawu[min_index])
        CDwu_t.append(CDwu[min_index])
        delta_array_wing_tail.append(delta_array_wing[min_index, :])

    CMwu_t = np.array(CMwu_t)
    CLwu_t = np.array(CLwu_t)
    CDwu_t = np.array(CDwu_t)
    alphawu_t = np.array(alphawu_t)
    delta_array_wing_tail = np.array(delta_array_wing_tail)

    delta_array_tail = get_streamline_tail(CMu, CMwu_t, CLwu_t, alphau, aircraft, tunnel)
    for i,point in enumerate(polar.points):
        point.et = esb + ewb[i] + ess[i]
        point.qInf = point.qInf * ((1 + point.et) ** 2)
        point.U = point.U * (1 + point.et)
        point.get_reynolds(aircraft)
        point.get_advanceratio(aircraft)
        alpha_temp = point.alpha
        point.alpha = point.alpha + delta_array_wing_tail[i, 0]
        point.alpha_t = alpha_temp + delta_array_tail[i, 0]
        point.CFL_w = CLwu_t[i] * (1 / (1 + point.et) ** 2)
        point.CFd_w = CDwu_t[i] * (1 / (1 + point.et) ** 2) + delta_array_wing_tail[i, 1]
        point.CMp_w = CMwu_t[i] * (1 / (1 + point.et) ** 2) + delta_array_wing_tail[i, 2]
        point.CFd = point.CFd * (1 / (1 + point.et) ** 2) + delta_array_wing_tail[i, 1]
        point.CMp = point.CMp * (1 / (1 + point.et) ** 2) + delta_array_wing_tail[i, 2] + delta_array_tail[i, 1]
        point.CFl = point.CFl * (1 / (1 + point.et) ** 2)
        point.CFl_t = point.CFl - point.CFL_w
        point.CFd_t = point.CFd - point.CFd_w
        point.CMp_t = point.CMp - point.CMp_w
        q_ratio = np.sqrt(get_slipstream(point.Tc,aircraft))
        delta_CL = point.CFl_t * ((1/q_ratio)**2 -1)
        point.CFl_t_0 = point.CFl_t + delta_CL

def interpolateTc(polar):
    J16_alpha = [0.0, 4.0]
    J16_Tc = [0.318, 0.324]
    J16_Tc_fit = np.interp(3, J16_alpha, J16_Tc)

    J20_alpha = [0.0, 4.0]
    J20_Tc = [0.1015, 0.1127]
    J20_Tc_fit = np.interp(3, J20_alpha, J20_Tc)

    J24_alpha = [0.0, 4.0]
    J24_Tc = [-0.0117, -0.0056]
    J24_Tc_fit = np.interp(3, J24_alpha, J24_Tc)

    Juc = np.array([1.6, 2.0, 2.4])
    Tc = np.array([J16_Tc_fit, J20_Tc_fit, J24_Tc_fit])

    for i, point in enumerate(polar.points):
        point.Tc = Tc[i]

def correct_deflection(polar,polar_tailOff,polar_J16u,polar_J20u,polar_J24u,polar_J16,polar_J20,polar_J24,aircraft,tunnel):

    'Uncorrected tail off'
    CLwu = np.array([point.CFl for point in polar_tailOff.points])
    alphawu = np.array([point.alpha for point in polar_tailOff.points])
    CMwu = np.array([point.CMp for point in polar_tailOff.points])
    CDwu = np.array([point.CFd for point in polar_tailOff.points])
    delta_array_wing = get_streamline_wing(CLwu, alphawu, aircraft, tunnel)
    delta_wing_alpha = np.interp(3.0,alphawu,delta_array_wing[:,0])
    delta_wing_drag = np.interp(3.0, alphawu, delta_array_wing[:,1])
    delta_wing_moment = np.interp(3.0, alphawu, delta_array_wing[:,2])
    CLwu_val = np.interp(3.0,alphawu,CLwu)
    CDwu_val = np.interp(3.0,alphawu,CDwu)
    CMwu_val = np.interp(3.0,alphawu,CMwu)

    for i, point in enumerate(polar.points):
        if i == 0:
            polar_tailOn = polar_J16u
            et_polar = polar_J16
        elif i == 1:
            polar_tailOn= polar_J20u
            et_polar = polar_J20
        elif i == 2:
            polar_tailOn = polar_J24u
            et_polar = polar_J24

        'Uncorrected net'
        CLu = np.array([point.CFl for point in polar_tailOn.points])
        CDu = np.array([point.CFd for point in polar_tailOn.points])
        CMu = np.array([point.CMp for point in polar_tailOn.points])
        alphau = np.array([point.alpha for point in polar_tailOn.points])

        'Get same alphas for uncorrected wing'
        CMwu_t = []
        CLwu_t = []
        alphawu_t = []
        CDwu_t = []
        delta_array_wing_tail = []
        for i, alpha in enumerate(alphau):
            min_index = np.argmin(np.abs(alphawu - alpha))
            CMwu_t.append(CMwu[min_index])
            CLwu_t.append(CLwu[min_index])
            alphawu_t.append(alphawu[min_index])
            CDwu_t.append(CDwu[min_index])
            delta_array_wing_tail.append(delta_array_wing[min_index, :])

        CMwu_t = np.array(CMwu_t)
        CLwu_t = np.array(CLwu_t)
        CDwu_t = np.array(CDwu_t)
        alphawu_t = np.array(alphawu_t)

        delta_array_tail = get_streamline_tail(CMu, CMwu_t, CLwu_t, alphau, aircraft, tunnel)
        delta_tail_moment = np.interp(3.0,alphau,delta_array_tail[:,1])

        et_alpha = np.array([point.alpha for point in et_polar.points])
        et_value = np.array([point.et for point in et_polar.points])
        point.et = np.interp(3.0,et_alpha,et_value)
        point.qInf = point.qInf * ((1 + point.et) ** 2)
        point.U = point.U * (1 + point.et)
        point.get_reynolds(aircraft)
        point.get_advanceratio(aircraft)
        alpha_temp = point.alpha
        point.alpha = point.alpha + delta_wing_alpha
        point.alpha_t = alpha_temp
        point.CFL_w = CLwu_val * (1 / (1 + point.et) ** 2)
        point.CFd_w = CDwu_val * (1 / (1 + point.et) ** 2) + delta_wing_drag
        point.CMp_w = CMwu_val * (1 / (1 + point.et) ** 2) + delta_wing_moment
        point.CFd = point.CFd * (1 / (1 + point.et) ** 2) + delta_wing_drag
        point.CMp = point.CMp * (1 / (1 + point.et) ** 2) + delta_wing_moment + delta_tail_moment
        point.CFl = point.CFl * (1 / (1 + point.et) ** 2)
        point.CFl_t = point.CFl - point.CFL_w
        point.CFd_t = point.CFd - point.CFd_w
        point.CMp_t = point.CMp - point.CMp_w
