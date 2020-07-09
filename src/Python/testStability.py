from pythonwrapper import shell
from pythonwrapper import impactDataIndex, angleDataIndex
from pythonwrapper import postPenDataIndex
import numpy as np
import matplotlib
matplotlib.use('TkAgg',warn=False, force=True)
import matplotlib.pyplot as plt



s = shell(.460, 780, .292, 1460, 2574, 6, .033, 76, 45, 60, "Yamato")
impacts = {}

s.setDtMin(.1)
#s.calcImpactForwardEuler()
#impacts['Forward Euler .1'] = s.getImpact()

#s.calcImpactRungeKutta2()
#impacts['RungeKutta2 .1'] = s.getImpact()

#s.calcImpactRungeKutta4()
#impacts['RungeKutta4 .1'] = s.getImpact()


s.setDtMin(.04)
#s.calcImpactForwardEuler()
#impacts['Forward Euler .04'] = s.getImpact()

#s.calcImpactRungeKutta2()
#impacts['RungeKutta2 .04'] = s.getImpact()

s.calcImpactRungeKutta4()
#impacts['RungeKutta4 .04'] = s.getImpact()

s.setDtMin(.02)
s.calcImpactAdamsBashforth5()
impacts['AdamsBashforth 5 .02'] = s.getImpact()

#s.calcImpactForwardEuler()
#impacts['Forward Euler .02'] = s.getImpact()

#s.calcImpactRungeKutta2()
#impacts['RungeKutta2 .02'] = s.getImpact()

#s.calcImpactRungeKutta4()
#impacts['RungeKutta4 .02'] = s.getImpact()

s.setDtMin(.01)
s.calcImpactAdamsBashforth5()
impacts['AdamsBashforth 5 .01'] = s.getImpact()

s.calcImpactForwardEuler()
impacts['Forward Euler .01'] = s.getImpact()

s.calcImpactRungeKutta4()
impacts['RungeKutta4 .01'] = s.getImpact()

#s.calcImpactRungeKutta4Hybrid()
#impacts['RungeKutta4Hybrid .1/.01'] = s.getImpact()

s.setDtMin(.001)
s.calcImpactRungeKutta4()
impacts['RungeKutta4 .001'] = s.getImpact()

s.calcImpactForwardEuler()
impacts['Forward Euler .001'] = s.getImpact()

s.calcImpactRungeKutta2()
impacts['RungeKutta2 .001'] = s.getImpact()

s.setDtMin(.0001)
s.calcImpactRungeKutta4()
reference = s.getImpact()

s.calcImpactForwardEuler()
impacts['Forward Euler .0001'] = s.getImpact()

s.calcImpactRungeKutta2()
impacts['RungeKutta2 .0001'] = s.getImpact()

ax1 = plt.subplot(211)
ax2 = plt.subplot(212)
#ax3 = plt.subplot(213)

ax1L = []
ax2L = []
#ax3L = []

for key, value in impacts.items():
    ax1V, = ax1.plot(
        value[impactDataIndex.launchA,:], 
        value[impactDataIndex.distance,:] - reference[impactDataIndex.distance,:],
        label=key)
    ax1L.append(ax1V)
    ax2V, = ax2.plot(
        value[impactDataIndex.launchA,:], 
        value[impactDataIndex.tToTarget,:] - reference[impactDataIndex.tToTarget,:],
        label=key)
    ax2L.append(ax2V)
    '''ax3V, = ax3.plot(
        value[impactDataIndex.launchA,:], 
        value[impactDataIndex.ePenHN,:] - reference[impactDataIndex.ePenHN,:],
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

