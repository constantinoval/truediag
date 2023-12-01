# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 10:37:45 2018

@author: Konstantinov
"""

import numpy as np
import matplotlib.pyplot as plt
import shutil
import os
from copy import deepcopy

prepost=r'C:\LSTC\LS-PrePost\4.2-win32\lsprepost4.2_win32.exe'
procRezCmd="""openc d3plot "{0}\\step-n\\d3plot"
ac
genselect target element
genselect target part
genselect part add part 1/0
maxvalue 7
xyplot 1 savefile xypair "ep" 1 all
switch2pick
ascii nodfor open "{0}\\step-n\\nodfor" 0
ascii nodfor changelist 2
ascii nodfor plot 7 1
xyplot 1 filter none  60.00 msec 9
xyplot 1 savefile xypair "f" 1 all
"""
f=open('procRez.cfile','w')
f.write(procRezCmd.format(os.path.abspath(os.curdir)))
f.close()

expData=np.loadtxt('exp_data.txt', skiprows=1, usecols=[0,1,2], unpack=True,
                   delimiter=',')
expData[2]*=-2*np.pi

expF=lambda t: np.interp(t, expData[0], expData[2])

f=open('model/bcs.k', 'w')
f.write('*keyword\n')
f.write('*define_curve\n1\n')
for t, v in zip(*expData[0:2]):
    f.write('%f, %f\n' % (t,v))
f.write('*control_termination\n%f\n' % (expData[0][-1]))
f.write('*end')
f.close()
tt=np.linspace(0.1,expData[0][-1], 100)


sig0=250
epmax=0.6

def integrate(y, x):
    rez=[0]
    for i in range(1,len(y)):
        rez.append(rez[-1]+0.5*(y[i]+y[i-1])*(x[i]-x[i-1]))
    return np.array(rez)

def calcDiagr(expData, L0=10, D0=5):
    de=expData[1]/L0
    e=integrate(de, expData[0])
    s=expData[2]/np.pi/D0**2*4
    s=s*(1+e)
    e=np.log(1+e)
    return e,s

def fitt(k, x):
    rez=[]
    n=len(k)
    for xx in x:
        rr=0
        for i, pp in enumerate(k):
            rr+=pp*xx**(n-i-1)
        rez.append(rr)
    return rez

#plt.plot(expData[0], expData[2])
e0, s0 = calcDiagr(expData)
for i, ss in enumerate(s0):
    if ss>=sig0:
        break
istart=i
for i, ee in enumerate(e0):
    if ee>=epmax:
        break
iend=i
e0=e0[istart:iend]
e0-=e0[0]
s0=s0[istart:iend]
einit=e0
sinit=s0
e0=[0]
s0=[250]

f=plt.figure()
ax=f.add_subplot(111)
ax.plot(e0, s0)
f2=plt.figure()
ax2=f2.add_subplot(111)
f3=plt.figure()
ax3=f3.add_subplot(111)
ax2.plot(expData[0], expData[2])
curve=[einit, sinit]
curvef=lambda e: np.interp(e, *curve)

for i in range(10):
    if os.path.exists('step-n'):
        shutil.rmtree('step-n')
    os.mkdir('step-n')
    shutil.copyfile('model/main.k', 'step-n/main.k')
    shutil.copyfile('model/bcs.k', 'step-n/bcs.k')
    shutil.copyfile('model/solid.k', 'step-n/solid.k')
    f=open('step-n/material.k', 'w')
    f.write("""*keyword
*mat_piecewise_linear_plasticity
1, 7.850e-3, 2.E+5, 0.28
0,0,100


*define_curve
100
""")
    for ee, ss in zip(*curve):
        f.write("{:e}, {:e}\n".format(ee,ss))
    f.write("{:e}, {:e}\n".format(100,curve[1][-1]))
    f.write('*end')
    f.close()
    os.system('cd step-n && run_main')
    os.system(prepost+' -nographics c=procRez.cfile')
    fn=np.loadtxt('f', skiprows=1, unpack=True)
    fn[1]*=-2*np.pi
    ep=np.loadtxt('ep', skiprows=1, unpack=True)
    ff=lambda t: np.interp(t, fn[0], fn[1])
    eep=lambda t: np.interp(t, ep[0], ep[1])
    beta=expF(tt)/ff(tt)
    ax3.plot(beta)
    enew1=list(e0)
    snew1=list(s0)
    for j in range(100):
        eee=eep(tt[j])
        if j==0 and not len(e0):
            enew1.append(eee)
            snew1.append(beta[j]*curvef(eee))
        if eee<=enew1[-1]:
            continue
        enew1.append(eee)
        snew1.append(beta[j]*curvef(eee))
    enew=deepcopy(enew1)
    snew=deepcopy(snew1)
    p=np.polyfit(enew, snew,5)
    snew=fitt(p, enew)
    curve=[enew, snew]
    curvef=lambda e: np.interp(e, enew, snew)

    ax.plot(enew,snew)
    ax2.plot(fn[0], fn[1])
    shutil.copyfile('f', 'f-step{}'.format(i))

    f=open('curve{}'.format(i), 'w')
    for ee, ss in zip(*curve):
        f.write("{:e}, {:e}\n".format(ee,ss))
    f.close()

plt.show()
