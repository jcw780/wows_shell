from pythonwrapper import shell
from pythonwrapper import impactDataIndex, angleDataIndex
from pythonwrapper import postPenDataIndex
import numpy as np
import matplotlib
matplotlib.use('TkAgg',warn=False, force=True)
import matplotlib.pyplot as plt



s = shell(.460, 780, .292, 1460, 2574, 6, .033, 76, 45, 60, "Yamato")

s.calcImpactForwardEuler()
fE = s.getImpact()

s.calcImpactRungeKutta()
rk = s.getImpact()

s.calcImpactRungeKuttaHybrid()
rkH = s.getImpact()

s.setDtMin(.04)
s.calcImpactRungeKutta()
rkF = s.getImpact()

s.setDtMin(.0001)
s.calcImpactRungeKutta()
rkF2 = s.getImpact()


plt.subplot(211)
d1, = plt.plot(
    rkF2[int(impactDataIndex.launchA), :], 
    rkH[int(impactDataIndex.distance), :] - rkF2[int(impactDataIndex.distance), :], 
    label="rkH")
d2, = plt.plot(
    rk[int(impactDataIndex.launchA), :], 
    rkF[int(impactDataIndex.distance), :] - rkF2[int(impactDataIndex.distance), :], 
    label="rkF")
d3, = plt.plot(
    rk[int(impactDataIndex.launchA), :], 
    fE[int(impactDataIndex.distance), :] - rkF2[int(impactDataIndex.distance), :], 
    label="FE")
d4, = plt.plot(
    rk[int(impactDataIndex.launchA), :], 
    rk[int(impactDataIndex.distance), :] - rkF2[int(impactDataIndex.distance), :], 
    label="rk")
plt.legend(handles=[d1, d2, d3, d4])
plt.subplot(212)
t1, = plt.plot(
    rk[int(impactDataIndex.launchA), :], 
    rkH[int(impactDataIndex.tToTarget), :] - rkF2[int(impactDataIndex.tToTarget), :], 
    label="rkH")
t2, = plt.plot(
    rk[int(impactDataIndex.launchA), :], 
    rkF[int(impactDataIndex.tToTarget), :] - rkF2[int(impactDataIndex.tToTarget), :], 
    label="rkF")
t3, = plt.plot(
    rk[int(impactDataIndex.launchA), :], 
    fE[int(impactDataIndex.tToTarget), :] - rkF2[int(impactDataIndex.tToTarget), :], 
    label="FE")
t4, = plt.plot(
    rk[int(impactDataIndex.launchA), :], 
    rk[int(impactDataIndex.tToTarget), :] - rkF2[int(impactDataIndex.tToTarget), :], 
    label="rk")
plt.legend(handles=[t1, t2, t3, t4])
plt.show()