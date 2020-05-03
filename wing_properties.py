import numpy as np

b = 1.397
S = 0.2172
cmac = 0.165
A = 8.98
taper = 0.40
dihedral = 4
twist = 2
t = 0.15
Q =0.0030229

tunnel_b = 1800/1000.
tunnel_h = 1250/1000.
tunnel_S = tunnel_b*tunnel_h - (0.5*300*300*4)/1000000.
tunnel_A = tunnel_b/tunnel_h

f1w = b/tunnel_b
t1w = 0.87
k1w = 1.02
esb_w = (k1w*t1w*Q)/((tunnel_S)**(3/2.))

l = 1.342
d = 0.142
f_S = np.pi*(d/2)**2
f_V = 0.0160632
t1f = 0.86
k3f = 0.915
esb_f = (k3f*t1f*f_V)/((tunnel_S)**(3/2.))

t1ss = 0.86
k1ss = 0.90
ss_V = 0.0035296
esb_ss = (k1ss*t1ss*ss_V)/((tunnel_S)**(3/2.))

t1ms = 0.86
k1ms = 0.90
ms_V = 0.0004491
esb_ms = (t1ms*k1ms*ms_V)/((tunnel_S)**(3/2.))




