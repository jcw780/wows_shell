from pythonwrapper import shell
from pythonwrapper import impactDataIndex, angleDataIndex
from pythonwrapper import postPenDataIndex
import numpy as np
import matplotlib
matplotlib.use('TkAgg',warn=False, force=True)
import matplotlib.pyplot as plt
import os
#os.environ['KMP_DUPLICATE_LIB_OK']='True'


s = shell(.460, 780, .292, 1460, 2574, 6, .033, 76, 45, 60, "Yamato")
s.calcImpact()
s.calcAngles(70, 0)
angles = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 ]
s.calcPostPen(70, 0, angles, True, False)

impact = s.getImpact()
angle = s.getAngles()
postPen = s.getPostPen()
print(postPen.shape)
s.printAngles()
print(angle)

'''
plt.plot(impact[int(impactDataIndex.distance),:], impact[int(impactDataIndex.rawPen),:])
plt.plot(impact[int(impactDataIndex.distance),:], impact[int(impactDataIndex.ePenH),:])
plt.plot(impact[int(impactDataIndex.distance),:], impact[int(impactDataIndex.ePenHN),:])

'''
'''

plt.plot(angle[int(angleDataIndex.distance), :], angle[int(angleDataIndex.fuseD), :])
plt.plot(angle[int(angleDataIndex.distance), :], angle[int(angleDataIndex.armorD), :])
plt.plot(angle[int(angleDataIndex.distance), :], angle[int(angleDataIndex.ra0D), :])
plt.plot(angle[int(angleDataIndex.distance), :], angle[int(angleDataIndex.ra1D), :])
'''

'''
for i in range(len(angles)):
    plt.plot(postPen[int(postPenDataIndex.distance),i,:], postPen[int(postPenDataIndex.xwf),i,:])
    plt.plot(postPen[int(postPenDataIndex.distance),i,:], postPen[int(postPenDataIndex.x),i,:])
'''
plt.show()

