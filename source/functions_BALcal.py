import numpy as np

def get_BALcalfactors():
    p = [0.009200961, 0.01594834, 0.06134184, 0.06143589, 0.02461131, 0.01231626]
    p = np.array(p)

    c = [[1, -0.4414104e-03, -0.8049510e-04, -0.1115902e-03, -0.1456587e-03, -0.6872960e-03],
        [0.1102029e-03, 1, -0.1138073e-03, 0.5521437e-04, 0, 0],
        [-0.2580602e-03, -0.3546250e-03, 1, 0, 0, 0],
        [-0.2580602e-03, -0.3546250e-03, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, -0.1122973e-02],
        [0, 0, 0, 0, 0.1084830e-03, 1]]
    c = np.array(c)

    arm = [0.001923, 0.000250, -0.000209, -0.525748, -0.525765, 0.000114, -0.262875, 0.263078, -0.000975, -1.050331,
           -1.050106, -1.049434]
    arm = np.array(arm)

    phiYbal = 0
    phiZbal = 0

    FX_cor = 0

    x_bend = 0
    y_bend = 0

    f = open('source/eij-bal','r')
    lines = f.readlines()
    f.close()
    e = []
    for line in lines:
        line = line.split()
        e.append(float(line[0]))

    e = np.array(e)
    e = e.reshape((12,6))

    return p,c,arm,FX_cor,x_bend,y_bend,e

def zero_BALdata(BAL_class):
    index = BAL_class.index
    index0 = BAL_class.index0
    data = BAL_class.data
    data0 = BAL_class.data0
    data_temp = np.zeros((len(data),6))
    for i,data_i in enumerate(data):
        alpha_i = data_i[index.AoA]
        dalpha_i = np.abs(data0[:,index0.AoA]-alpha_i)
        min_index0 = np.argmin(dalpha_i)
        data_temp[i,:] = data_i[index.B1:index.B6+1] - data0[min_index0,index0.B1:index0.B6+1]
    return data_temp


def process_BALcallibration(BAL_class):
    p, c, arm, FX_cor, x_bend, y_bend, e = get_BALcalfactors()
    data_temp = zero_BALdata(BAL_class)
    index = BAL_class.index
    F1 = []
    F2 = []
    F3 = []
    M1 = []
    M2 = []
    M3 = []
    B1_cal = []
    B2_cal = []
    B3_cal = []
    B4_cal = []
    B5_cal = []
    B6_cal = []
    for data_i in data_temp:
        bal_vec = data_i[:]
        f = np.matmul(c,(p*bal_vec))
        B1_cal.append(bal_vec[0])
        B2_cal.append(bal_vec[1])
        B3_cal.append(bal_vec[2])
        B4_cal.append(bal_vec[3])
        B5_cal.append(bal_vec[4])
        B6_cal.append(bal_vec[5])
        F1.append(f[0])
        F2.append(f[1]+f[5])
        F3.append(f[2]+f[3]+f[4])

        arm_new = arm + np.matmul(e,f)

        M1.append(-f[1]*arm_new[10]-f[5]*arm_new[11]+f[2]*arm_new[6]+f[3]*arm_new[7]+f[4]*arm_new[8])
        M2.append(f[0]*arm_new[9]-f[2]*arm_new[1]-f[3]*arm_new[2]-f[4]*arm_new[3])
        M3.append(-f[0]*arm_new[5]+f[1]*arm_new[0]+f[5]*arm_new[4])

    B1_cal = np.array(B1_cal).reshape((len(B1_cal),1))
    B2_cal = np.array(B2_cal).reshape((len(B2_cal),1))
    B3_cal = np.array(B3_cal).reshape((len(B3_cal),1))
    B4_cal = np.array(B4_cal).reshape((len(B4_cal),1))
    B5_cal = np.array(B5_cal).reshape((len(B5_cal),1))
    B6_cal = np.array(B6_cal).reshape((len(B6_cal),1))
    F1 = np.array(F1).reshape((len(F1),1))
    F2 = np.array(F2).reshape((len(F2),1))
    F3 = np.array(F3).reshape((len(F3),1))
    M1 = np.array(M1).reshape((len(M1),1))
    M2 = np.array(M2).reshape((len(M2),1))
    M3 = np.array(M3).reshape((len(M3),1))
    data_cal = np.concatenate((B1_cal, B2_cal, B3_cal, B4_cal, B5_cal, B6_cal, F1, F2, F3, M1, M2, M3), axis=1)
    return data_cal