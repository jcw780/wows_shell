from pythonwrapper import shell
from pythonwrapper import impactIndices, angleIndices, postPenIndices
import numpy as np
import matplotlib
matplotlib.use('TkAgg',warn=False, force=True)
import matplotlib.pyplot as plt
import os
#os.environ['KMP_DUPLICATE_LIB_OK']='True'


s = shell(.460, 780, .292, 1460, 2574, 6, .033, 76, 45, 60, 0, "Yamato")
s.calcImpactRungeKutta4()
s.calcAngles(70, 0)
angles = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 ]
s.calcPostPen(70, 0, angles, True, False)

impact = s.getImpact()
angle = s.getAngles()
postPen = s.getPostPen()
print(postPen.shape)
s.printAngles()
print(angle)

plt.subplot(311)
plt.plot(impact[int(impactIndices.distance),:], impact[int(impactIndices.rawPenetration),:])
plt.plot(impact[int(impactIndices.distance),:], impact[int(impactIndices.effectivePenetrationHorizontal),:])
plt.plot(impact[int(impactIndices.distance),:], impact[int(impactIndices.effectivePenetrationHorizontalNormalized),:])

plt.subplot(312)
plt.plot(angle[int(angleIndices.distance), :], angle[int(angleIndices.fuseDegrees), :])
plt.plot(angle[int(angleIndices.distance), :], angle[int(angleIndices.armorDegrees), :])
plt.plot(angle[int(angleIndices.distance), :], angle[int(angleIndices.ricochetAngle0Degrees), :])
plt.plot(angle[int(angleIndices.distance), :], angle[int(angleIndices.ricochetAngle1Degrees), :])

plt.subplot(313)
for i in range(len(angles)):
    plt.plot(postPen[int(postPenIndices.distance),i,:], postPen[int(postPenIndices.xwf),i,:])
    plt.plot(postPen[int(postPenIndices.distance),i,:], postPen[int(postPenIndices.x),i,:])

plt.show()

