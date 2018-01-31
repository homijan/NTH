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
Te = 10000.0
gradTe = -1.0
Zbar = 47.0
sigma = 1e15
# The heat flux after integration takes the form
# qH = me/Zbar/sigma*128/(2*pi)**0.5*(kB/me)**(7/2)*T**(5/2)*gradT,
# where mfp_ei = v**4/sigma/Zbar, i.e. sigma corresponds to ee collisions.
# Physical fix by a magic constant.
# Z = 1
#cmag = 1.882
# Z = 5
#cmag = 1.59575
# Z = 10
#cmag = 1.565
# Z = 20
#cmag = 1.56
# Z = 50
#cmag = 1.585
# Z = 100
#cmag = 1.65
# Z = 200
#cmag = 1.8
#coeffs = np.array([0.5313496280552604, 0.6266645777847407, 0.6389776357827476, 0.641025641025641, 0.6309148264984227, 0.6060606060606061, 0.5555555555555556])
#Zbars = [1.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0]
#p1 = 688.9
#p2 = 114.4
#q1 = 1038
#q2 = 474.1
#icoeffs = np.array([(p1*Zbars[i] + p2)/(Zbars[i]**2.0 + q1*Zbars[i] + q2) for i in range(7)])
#print icoeffs
#Zbar = Zbars[0]
#corr = coeffs[0]
#Zbar = Zbars[1]
#corr = coeffs[1]
#Zbar = Zbars[2]
#corr = coeffs[2]
#Zbar = Zbars[3]
#corr = coeffs[3]
#Zbar = Zbars[4]
#corr = coeffs[4]
#Zbar = Zbars[5]
#corr = coeffs[5]
#Zbar = Zbars[6]
#corr = coeffs[6]
corr = (688.9*Zbar + 114.4)/(Zbar**2.0 + 1038*Zbar + 474.1)
print("Zbar, corr:", Zbar, corr)
cmag = 1./corr
#Efield = 0.0
Efield = vTh(Te)**2.0*2.5*gradTe/Te

Nexpl = 5000
vexpl = np.linspace(ml_max*vTh(Te), ml_min*vTh(Te), Nexpl)
dvexpl = (ml_max - ml_min)*vTh(Te)/(Nexpl-1)

Nimpl = 5000
vimpl = np.linspace(ml_max*vTh(Te), ml_min*vTh(Te), Nimpl)
dvimpl = (ml_max - ml_min)*vTh(Te)/(Nimpl-1)

sol_expl = solve_fweuler(vexpl, 0.0, Te, gradTe, Zbar, Efield)
sol_impl = solve_bweuler(vimpl, 0.0, Te, gradTe, cmag*Zbar, Efield)

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
    # The mean free path has standard v^4 dependence, sigma is cross section
    # given by model and Zbar increases the effect of Coulomb potential in
    # ei collisions
    mfp = vp**4.0/sigma/Zbar
    dv = dvimpl
    #BGKf1[i] = - Zbar/(Zbar + 1.0)*((vp**2.0/2.0/vTh(Te)**2.0 - 1.5)*gradTe/Te - Efield/vTh(Te)**2.0)*fM(vp, Te)*vp*vp
    BGKf1[i] = - (Zbar + 0.24)/(Zbar + 4.2)*((vp**2.0/2.0/vTh(Te)**2.0 - 1.5)*gradTe/Te - Efield/vTh(Te)**2.0)*fM(vp, Te)*vp*vp
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
print("SHJ: ", BGKJ)
print("SHQ: ", BGKQ)
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
#plt.plot(Zbars, coeffs)
#plt.plot(Zbars, icoeffs)
#plt.plot(vexpl, M1q_expl, 'g--', label='M1q_expl')
plt.plot(vimpl/vTh(Te), M1q_impl, 'r--', label='M1q_impl')
plt.plot(vimpl/vTh(Te), BGKq, 'b', label='SHq')
#plt.plot(v, nM1q, 'b', label='nM1q')
#plt.plot(v, nBGKq, 'g', label='nBGKq')
#plt.plot(v, nM1j, 'b--', label='nM1j')
#plt.plot(v, nBGKj, 'g--', label='nBGKj')
plt.legend(loc='best')
plt.xlabel('v/vT')
plt.grid()
plt.show()
