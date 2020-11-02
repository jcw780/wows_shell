from pythonwrapper import shell, shellCalc
from pythonwrapper import impactIndices, angleIndices, postPenIndices
import numpy as np
import matplotlib
#matplotlib.use('TkAgg',warn=False, force=True)
import matplotlib.pyplot as plt

s = shell(.460, 780, .292, 1460, 2574, 6, .033, 76, 45, 60, 0, "Yamato")
impacts = {}

c = shellCalc()
c.setDtMin(.01)
c.calcImpactForwardEuler(s)
impacts['Forward Euler .01'] = s.getImpact()

c.setDtMin(.04)
c.calcImpactRungeKutta2(s)
impacts['Runge Kutta .01'] = s.getImpact()

c.setDtMin(.0001)
c.calcImpactRungeKutta4(s)
reference = s.getImpact()

c.calcImpactForwardEuler(s)
impacts['Forward Euler .0001'] = s.getImpact()

c.calcImpactRungeKutta2(s)
impacts['RungeKutta2 .0001'] = s.getImpact()

plt.plot(impacts['Forward Euler .01'][impactIndices.distance,:], impacts['Forward Euler .01'][impactIndices.effectivePenetrationHorizontalNormalized,:])
plt.plot(impacts['Runge Kutta .01'][impactIndices.distance,:], impacts['Runge Kutta .01'][impactIndices.effectivePenetrationHorizontalNormalized,:])
plt.plot(reference[impactIndices.distance,:], reference[impactIndices.effectivePenetrationHorizontalNormalized,:])

plt.show()