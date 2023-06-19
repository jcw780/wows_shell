from wows_shell import *
import numpy as np
import matplotlib
# matplotlib.use('TkAgg',warn=False, force=True)
import matplotlib.pyplot as plt
import os
# os.environ['KMP_DUPLICATE_LIB_OK']='True'
caliber = 0.457
muzzle = 732
drag = 0.37
mass = 1746
krupp = 2400
normal = 6
name = "Vermont"

s = shell(shellParams(caliber, muzzle, drag, mass,
                      krupp, normal, .033, 76, 45, 60, 0),
          dispersionParams(10, 2.8, 1000, 5000, 0.5, 0.2, 0.6, 0.8, 26630, 2.1), name)
c = shellCalc()
c.calcImpactForwardEuler(s)

impact = s.getImpact()
velocity = impact[int(impactIndices.impactVelocity), :]

updatedPenetration = krupp * (((velocity ) ** 2 * mass) ** 0.69) * (caliber ** -1.07) * 0.0000001

plt.title(f'{name} Regression vs. RE')
plt.plot(impact[int(impactIndices.distance), :],
         impact[int(impactIndices.rawPenetration), :], label='Regression')

plt.plot(impact[int(impactIndices.distance), :],
         updatedPenetration, label='RE')
plt.legend()
plt.show()