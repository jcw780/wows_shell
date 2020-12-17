from pythonwrapper import *

dataList = [
    (shellParams(.460, 780, .292, 1460, 2574,
                 6, .033, 76, 45, 60, 0),
     dispersionParams(10, 2.8, 1000, 5000,
                      0.5, 0.2, 0.6, 0.8,
                      26630, 2.1),
     "Yamato"),
    (shellParams(
        {"caliber": .457, "v0": 762, "cD": 0.2897,
         "mass": 1506, "krupp": 2485, "normalization": 6,
         "fuseTime": .033, "threshold": 76, "ricochet0": 45,
         "ricochet1": 60, "nonAP": 0}),
     dispersionParams(
         {"idealRadius": 10, "minRadius": 1.6, "idealDistance": 1000,
             "taperDistance": 5000,
             "delim": 0.5, "zeroRadius": 0.2, "delimRadius": 0.5, "maxRadius": 0.6,
             "maxDistance": 24250, "sigma": 1.8}),
     "Thunderer"),
    (shellParams(
        caliber=.406, v0=870, cD=0.2,
        mass=890, krupp=2890, normalization=6,
        fuseTime=.033, threshold=68, ricochet0=45,
        ricochet1=60, nonAP=0),
        dispersionParams(
            {"idealRadius": 8.5, "minRadius": 3.5, "idealDistance": 1000,
             "taperDistance": 5000,
             "delim": 0.3, "zeroRadius": 0.1, "delimRadius": 0.25, "maxRadius": 0.45,
             "maxDistance": 24680, "sigma": 1.9}),
        "Slava")
]

'''shellList = [
    shell(
        shellParams(.460, 780, .292, 1460, 2574,
                    6, .033, 76, 45, 60, 0),
        dispersionParams(10, 2.8, 1000, 5000,
                         0.5, 0.2, 0.6, 0.8,
                         26630, 2.1),
        "Yamato"
    ),
    shell(
        shellParams(
            {"caliber": .457, "v0": 762, "cD": 0.2897,
             "mass": 1506, "krupp": 2485, "normalization": 6,
             "fuseTime": .033, "threshold": 76, "ricochet0": 45,
             "ricochet1": 60, "nonAP": 0}),
        dispersionParams(
            {"idealRadius": 10, "minRadius": 1.6, "idealDistance": 1000,
             "taperDistance": 5000,
             "delim": 0.5, "zeroRadius": 0.2, "delimRadius": 0.5, "maxRadius": 0.6,
             "maxDistance": 24250, "sigma": 1.8}),
        "Thunderer"
    ),
    shell(
        shellParams(
            caliber=.406, v0=870, cD=0.2,
            mass=890, krupp=2890, normalization=6,
            fuseTime=.033, threshold=68, ricochet0=45,
            ricochet1=60, nonAP=0),
        dispersionParams(
            {"idealRadius": 8.5, "minRadius": 3.5, "idealDistance": 1000,
             "taperDistance": 5000,
             "delim": 0.3, "zeroRadius": 0.1, "delimRadius": 0.25, "maxRadius": 0.45,
             "maxDistance": 24680, "sigma": 1.9}),
        "Slava"
    )
]'''

shellList = []
for data in dataList:
    shellList.append(shell(*data))

c = shellCalc()
angles = [0, 10, 20]

for i, s in enumerate(shellList):
    s.setValues(*dataList[i])
    c.calcImpactForwardEuler(s)
    c.calcAngles(s, 70, 0)
    c.calcDispersion(s)
    c.calcPostPen(s, 70, 0, angles, True, False)

for s in shellList:
    # s.printImpact()
    # s.printAngles()
    # s.printDispersion()
    # s.printPostPen()
    print(s.getImpact(owned=False))
    print(s.getAngles(owned=False))
    print(s.getDispersion(owned=False))
    print(s.getPostPen(owned=False))
    print(generateShellHash(s))
