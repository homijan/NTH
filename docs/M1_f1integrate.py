import numpy as np
from scipy.integrate import odeint
from math import pi
from math import exp

kB = 1.0
me = 1.0
Te = 1.0

ne = 1.0

def vTh(T): 
    return (kB*T/me)**0.5

def fM(v, T):
    return ne/(vTh(T)**3.0*(2.0*pi)**1.5)*exp(-v**2.0/2.0/vTh(T)**2.0)

def rhs(v, T, gradT, Z, E):
    return Z/v*((v**2.0/2.0/vTh(T)**2.0 - 1.5)*gradT/T - E/vTh(T)**2.0)*fM(v, T)

def func(y, v, T, gradT, Z, E):
    f1, = y
    dydt = [(Z - 4.0)/v*f1 + rhs(v, T, gradT, Z, E)]
    return dydt

def solve_fweuler(v, f0, T, gradT, Z, E):
    N = len(v)
    f1 = np.zeros(N)
    f1[0] = f0
    for i in range(N-1):
        dv = v[i+1] - v[i]
        vp = v[i]
        #f1[i+1] = f1[i] + dv*rhs(vp, T, gradT, Z, E)
        f1[i+1] = (1.0 + dv*(Z - 4.0)/vp)*f1[i] + dv*rhs(vp, T, gradT, Z, E)
    return f1

def solve_bweuler(v, f0, T, gradT, Z, E):
    N = len(v)
    f1 = np.zeros(N)
    f1[0] = f0
    for i in range(N-1):
        dv = v[i+1] - v[i]
        vp = v[i]
        f1[i+1] = (f1[i] + dv*rhs(vp, T, gradT, Z, E))/(1.0 + dv*(4.0 - Z)/vp)
    return f1

#v0 = 6.4
#T0 = 92.0
#gT0 = 18.0
#dT = 1e-5
#print((v0**2.0/2.0/vTh(T0)**2.0-1.5)*fM(v0, T0)/T0)
#print((fM(v0, T0+dT)-fM(v0, T0))/dT)

# Optimal implicit solve
#N = 150
#ml = 7.0
#Te = 1000.0
#gradTe = -1.0
#Zbar = 1000.0

ml_max = 10.0
ml_min = 0.0
#ml_max = 6.5
#ml_min = 0.5
Te = 1000.0
gradTe = -1.0
Zbar = 100.0
#Efield = 0.0
Efield = vTh(Te)**2.0*2.5*gradTe/Te

Nexpl = 1000
vexpl = np.linspace(ml_max*vTh(Te), ml_min*vTh(Te), Nexpl)
dvexpl = (ml_max - ml_min)*vTh(Te)/(Nexpl-1)

Nimpl = 1000
vimpl = np.linspace(ml_max*vTh(Te), ml_min*vTh(Te), Nimpl)
dvimpl = (ml_max - ml_min)*vTh(Te)/(Nimpl-1)

sol_expl = solve_fweuler(vexpl, 0.0, Te, gradTe, Zbar, Efield)
sol_impl = solve_bweuler(vimpl, 0.0, Te, gradTe, Zbar, Efield)

# Post-process transport values
BGKf1 = np.zeros(Nimpl)
BGKj = np.zeros(Nimpl)
BGKq = np.zeros(Nimpl)
BGKJ = 0.0
BGKQ = 0.0
M1f1_impl = np.zeros(Nimpl)
M1j_impl = np.zeros(Nimpl)
M1q_impl = np.zeros(Nimpl)
M1J_impl = 0.0
M1Q_impl = 0.0
for i in range(Nimpl):
    vp = vimpl[i]
    mfp = vp**4.0
    dv = dvimpl
    BGKf1[i] = - Zbar/(Zbar + 1.0)*((vp**2.0/2.0/vTh(Te)**2.0 - 1.5)*gradTe/Te - Efield/vTh(Te)**2.0)*fM(vp, Te)*vp*vp
    BGKj[i] = mfp*vp*BGKf1[i]
    BGKq[i] = mfp*me/2.0*vp*vp*vp*BGKf1[i]
    BGKJ = BGKJ + 4.0*pi/3.0*BGKj[i]*dv
    BGKQ = BGKQ + 4.0*pi/3.0*BGKq[i]*dv
    M1f1_impl[i] = sol_impl[i]*vp*vp
    M1j_impl[i] = mfp*vp*M1f1_impl[i]
    M1q_impl[i] = mfp*vp*vp*vp*me/2.0*M1f1_impl[i]
    M1J_impl = M1J_impl + 4.0*pi/3.0*M1j_impl[i]*dv
    M1Q_impl = M1Q_impl + 4.0*pi/3.0*M1q_impl[i]*dv
M1f1_expl = np.zeros(Nexpl)
M1j_expl = np.zeros(Nexpl)
M1q_expl = np.zeros(Nexpl)
M1J_expl = 0.0
M1Q_expl = 0.0
for i in range(Nexpl):
    vp = vexpl[i]
    mfp = vp**4.0
    dv = dvexpl
    M1f1_expl[i] = sol_expl[i]*vp*vp
    M1j_expl[i] = mfp*vp*M1f1_expl[i]
    M1q_expl[i] = mfp*vp*vp*vp*me/2.0*M1f1_expl[i]
    M1J_expl = M1J_expl + 4.0*pi/3.0*M1j_expl[i]*dv
    M1Q_expl = M1Q_expl + 4.0*pi/3.0*M1q_expl[i]*dv

# Print integrated values
print("BGKJ: ", BGKJ)
print("BGKQ: ", BGKQ)
print("M1J_impl: ", M1J_impl)
print("M1Q_impl: ", M1Q_impl)
print("M1J_expl: ", M1J_expl)
print("M1Q_expl: ", M1Q_expl)
# Normalization for graphical output
#nBGKj = BGKj/max(abs(BGKj))
#nBGKq = BGKq/max(abs(BGKq))
#nM1j = M1j_impl/max(abs(M1j_impl))
#nM1q = M1q_impl/max(abs(M1q_impl))

import matplotlib.pyplot as plt
plt.plot(vexpl, M1q_expl, 'b', label='M1q_expl')
plt.plot(vimpl, M1q_impl, 'k--', label='M1q_impl')
plt.plot(vimpl, BGKq, 'g--', label='BGKq')
#plt.plot(v, nM1q, 'b', label='nM1q')
#plt.plot(v, nBGKq, 'g', label='nBGKq')
#plt.plot(v, nM1j, 'b--', label='nM1j')
#plt.plot(v, nBGKj, 'g--', label='nBGKj')
plt.legend(loc='best')
plt.xlabel('v')
plt.grid()
plt.show()
