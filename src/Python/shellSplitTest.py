from pythonwrapper import shellS, shellCalc
from pythonwrapper import impactIndices, angleIndices, postPenIndices

shellList = [
    shellS(.460, 780, .292, 1460, 2574, 6, .033, 76, 45, 60, 0, "Yamato"),
    shellS(.457, 762, .2897, 1506, 2485, 6, .033, 76, 45, 60.0, 0, "Thunderer")
]
c = shellCalc()
angles = [0, 10, 20, 30, 40, 50]

for s in shellList:
    c.calcImpactForwardEuler(s)
    c.calcAngles(s, 70, 0)
    c.calcPostPen(s, 70, 0, angles, True, False)

for s in shellList:
    s.printImpact()
    s.printAngles()
    s.printPostPen()
    print(s.getImpact())