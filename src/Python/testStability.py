from pythonwrapper import shell, shellCalc
from pythonwrapper import impactIndices, angleIndices, postPenIndices
import numpy as np
import matplotlib
#matplotlib.use('TkAgg',warn=False, force=True)
import matplotlib.pyplot as plt

#Yes, this can be written more efficiently

s = shell(.460, 780, .292, 1460, 2574, 6, .033, 76, 45, 60, 0, 'Yamato')
impacts = {}

c = shellCalc()
c.setDtMin(.1)
c.calcImpactForwardEuler(s)
impacts['Forward Euler .1'] = s.getImpact()

c.calcImpactRungeKutta2(s)
impacts['RungeKutta2 .1'] = s.getImpact()

c.calcImpactRungeKutta4(s)
impacts['RungeKutta4 .1'] = s.getImpact()


c.setDtMin(.04)
c.calcImpactForwardEuler(s)
impacts['Forward Euler .04'] = s.getImpact()

c.calcImpactRungeKutta2(s)
impacts['RungeKutta2 .04'] = s.getImpact()

c.calcImpactRungeKutta4(s)
impacts['RungeKutta4 .04'] = s.getImpact()

c.setDtMin(.02)
c.calcImpactAdamsBashforth5(s)
impacts['AdamsBashforth 5 .02'] = s.getImpact()

c.calcImpactForwardEuler(s)
impacts['Forward Euler .02'] = s.getImpact()

c.calcImpactRungeKutta2(s)
impacts['RungeKutta2 .02'] = s.getImpact()

c.calcImpactRungeKutta4(s)
impacts['RungeKutta4 .02'] = s.getImpact()

c.setDtMin(.01)
c.calcImpactAdamsBashforth5(s)
impacts['AdamsBashforth 5 .01'] = s.getImpact()

c.calcImpactForwardEuler(s)
impacts['Forward Euler .01'] = s.getImpact()

c.calcImpactRungeKutta4(s)
impacts['RungeKutta4 .01'] = s.getImpact()

c.setDtMin(.001)
c.calcImpactRungeKutta4(s)
impacts['RungeKutta4 .001'] = s.getImpact()

c.calcImpactForwardEuler(s)
impacts['Forward Euler .001'] = s.getImpact()

c.calcImpactRungeKutta2(s)
impacts['RungeKutta2 .001'] = s.getImpact()

c.setDtMin(.0001)
c.calcImpactRungeKutta4(s)
reference = s.getImpact()

c.calcImpactForwardEuler(s)
impacts['Forward Euler .0001'] = s.getImpact()

c.calcImpactRungeKutta2(s)
impacts['RungeKutta2 .0001'] = s.getImpact()

ax1 = plt.subplot(211)
ax2 = plt.subplot(212)
#ax3 = plt.subplot(213)

ax1L = []
ax2L = []
#ax3L = []

for key, value in impacts.items():
    ax1V, = ax1.plot(
        value[impactIndices.launchAngle,:], 
        value[impactIndices.distance,:] - reference[impactIndices.distance,:],
        label=key)
    ax1L.append(ax1V)
    ax2V, = ax2.plot(
        value[impactIndices.launchAngle,:], 
        value[impactIndices.timeToTarget,:] - reference[impactIndices.timeToTarget,:],
        label=key)
    ax2L.append(ax2V)
    '''ax3V, = ax3.plot(
        value[impactIndices.launchA,:], 
        value[impactIndices.ePenHN,:] - reference[impactIndices.ePenHN,:],
        label=key)
    ax3L.append(ax3V)'''

ax1.title.set_text('Distance error - rk4 .0001 reference')
ax1.set_xlabel('Launch Angle')
ax1.set_ylabel('Distance error (m)')
ax1.legend(handles=ax1L)
ax2.title.set_text('Flight time error - rk4 .0001 reference')
ax2.set_xlabel('Launch Angle')
ax2.set_ylabel('Time error (s)')
ax2.legend(handles=ax2L)

plt.show()

