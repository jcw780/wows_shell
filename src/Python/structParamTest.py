from pythonwrapper import *

shellList = [
    shell(
        shellParams(.460, 780, .292, 1460, 2574,
                    6, .033, 76, 45, 60, 0),
        dispersionParams(10, 2.8, 1000, 5000,
                         0.5, 0.2, 0.6, 0.8,
                         26630, 2.1),
        "Yamato"
    ),
    shell(
        shellParams(.457, 762, .2897, 1506, 2485,
                    6, .033, 76, 45, 60.0, 0),
        dispersionParams(10, 1.6, 1000, 5000,
                         0.5, 0.2, 0.5, 0.6,
                         24250, 1.8),
        "Thunderer"
    )
]
c = shellCalc()
angles = [0, 10, 20, 30, 40, 50]

for s in shellList:
    c.calcImpactForwardEuler(s)
    c.calcAngles(s, 70, 0)
    c.calcDispersion(s)
    c.calcPostPen(s, 70, 0, angles, True, False)

for s in shellList:
    s.printImpact()
    s.printAngles()
    s.printDispersion()
    s.printPostPen()
    print(s.getImpact())
    print(s.getAngles())
    print(s.getDispersion())
    print(s.getPostPen())
    print(generateShellHash(s))
