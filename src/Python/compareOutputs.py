from pythonwrapper import shell
from pythonwrapper import impactIndices, angleIndices
from pythonwrapper import postPenIndices
import numpy as np
import matplotlib
matplotlib.use('TkAgg',warn=False, force=True)
import matplotlib.pyplot as plt

s = shell(.460, 780, .292, 1460, 2574, 6, .033, 76, 45, 60, "Yamato")
impacts = {}

s.setDtMin(.01)
s.calcImpactForwardEuler()
impacts['Forward Euler .01'] = s.getImpact()

s.setDtMin(.04)
s.calcImpactRungeKutta2()
impacts['Runge Kutta .01'] = s.getImpact()

s.setDtMin(.0001)
s.calcImpactRungeKutta4()
reference = s.getImpact()

#s.calcImpactForwardEuler()
#impacts['Forward Euler .0001'] = s.getImpact()

#s.calcImpactRungeKutta2()
#impacts['RungeKutta2 .0001'] = s.getImpact()

plt.plot(impacts['Forward Euler .01'][impactIndices.distance,:], impacts['Forward Euler .01'][impactIndices.ePenHN,:])
plt.plot(impacts['Runge Kutta .01'][impactIndices.distance,:], impacts['Runge Kutta .01'][impactIndices.ePenHN,:])
plt.plot(reference[impactIndices.distance,:], reference[impactIndices.ePenHN,:])

plt.show()