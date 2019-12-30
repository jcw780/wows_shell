from pythonwrapper import shell
from pythonwrapper import impactDataIndex
from pythonwrapper import postPenDataIndex
import numpy as np
import matplotlib
matplotlib.use('TkAgg',warn=False, force=True)
import matplotlib.pyplot as plt
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'


s = shell(780, .460, 2574, 1460, 6, .292, "Yamato", 76.0, .033 )
s.calcImpact()
angles = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 ]
s.calcPostPen(70, angles)

impact = s.getImpact()
postPen = s.getPostPen()
print(postPen.shape)

'''
plt.plot(impact[int(impactDataIndex.distance),:], impact[int(impactDataIndex.rawPen),:])
plt.plot(impact[int(impactDataIndex.distance),:], impact[int(impactDataIndex.ePenH),:])
plt.plot(impact[int(impactDataIndex.distance),:], impact[int(impactDataIndex.ePenHN),:])
'''
for i in range(len(angles)):
    plt.plot(postPen[int(postPenDataIndex.distance),i,:], postPen[int(postPenDataIndex.xwf),i,:])
    plt.plot(postPen[int(postPenDataIndex.distance),i,:], postPen[int(postPenDataIndex.x),i,:])

plt.show()

