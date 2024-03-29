from wows_shell import *
import matplotlib.pyplot as plt

s = shell(shellParams(.460, 780, .292, 1460,
                      2574, 6, .033, 76, 45, 60, 0), "Yamato")

c = shellCalc()
c.setDtMin(.01)
c.setMax(90)

c.calcImpactForwardEuler(s, addTraj=True)

impactData = s.getImpact()
print(impactData.shape)

traj = s.getTrajectory(400)
print(traj)

plt.plot(traj[0,:], traj[1,:])
plt.plot(traj[0,:], traj[2,:])
plt.show()