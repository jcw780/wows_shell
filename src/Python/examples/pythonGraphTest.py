from wows_shell import *
import numpy as np
import matplotlib
# matplotlib.use('TkAgg',warn=False, force=True)
import matplotlib.pyplot as plt
import os
# os.environ['KMP_DUPLICATE_LIB_OK']='True'


s = shell(shellParams(.460, 780, .292, 1460,
                      2574, 6, .033, 76, 45, 60, 0),
          dispersionParams(10, 2.8, 1000, 5000, 0.5, 0.2, 0.6, 0.8, 26630, 2.1), "Yamato")
c = shellCalc()
c.calcImpactRungeKutta4(s)
c.calcAngles(s, 70, 0)
c.calcDispersion(s, int(verticalTypes.horizontal))
angles = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
c.calcPostPen(s, 70, 0, angles, True, False)

impact = s.getImpact()
angle = s.getAngles()
postPen = s.getPostPen()
dispersion = s.getDispersion()
print(postPen.shape)

ax1 = plt.subplot(411)
ax1.plot(impact[int(impactIndices.distance), :],
         impact[int(impactIndices.rawPenetration), :])
ax1.plot(impact[int(impactIndices.distance), :],
         impact[int(impactIndices.effectivePenetrationHorizontal), :])
ax1.plot(impact[int(impactIndices.distance), :], impact[int(
    impactIndices.effectivePenetrationHorizontalNormalized), :])

ax1.grid(b=True)

ax2 = ax1.twinx()
ax2.plot(impact[int(impactIndices.distance), :],
         impact[int(impactIndices.impactAngleHorizontalDegrees), :])

ax2.grid(b=True)

ax1 = plt.subplot(412)
ax1.plot(impact[int(impactIndices.distance), :],
         angle[int(angleIndices.fuseDegrees), :])
ax1.plot(impact[int(impactIndices.distance), :],
         angle[int(angleIndices.armorDegrees), :])
ax1.plot(impact[int(impactIndices.distance), :],
         angle[int(angleIndices.ricochetAngle0Degrees), :])
ax1.plot(impact[int(impactIndices.distance), :],
         angle[int(angleIndices.ricochetAngle1Degrees), :])
plt.grid(b=True)

ax1 = plt.subplot(413)
ax1.plot(impact[int(impactIndices.distance), :],
         dispersion[int(dispersionIndices.maxVertical), :])
ax1.plot(impact[int(impactIndices.distance), :],
         dispersion[int(dispersionIndices.standardVertical), :])
ax1.plot(impact[int(impactIndices.distance), :],
         dispersion[int(dispersionIndices.halfVertical), :])
plt.grid(b=True)

ax1 = plt.subplot(414)
for i in range(len(angles)):
    ax1.plot(impact[int(impactIndices.distance), :],
             postPen[int(postPenIndices.xwf), i, :])
    ax1.plot(impact[int(impactIndices.distance), :],
             postPen[int(postPenIndices.x), i, :])

plt.grid(b=True)
plt.show()
