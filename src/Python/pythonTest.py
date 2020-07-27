from pythonwrapper import shell, impactIndices, angleIndices, postPenIndices
import numpy as np
s = shell(.460, 780, .292, 1460, 2574, 6, .033, 76, 45, 60, 0, 'Yamato')
s.calcImpactForwardEuler()
print("Standard Done")
s.calcPostPen(400.0, 0, [0, 10, 20], True, False)
print("Postpen done")

n1 = s.getImpact()
print(n1.shape)

n1r = np.round(n1, 2)

for y in range(250):
    pS = F''
    for x in range(13):
        pS = F'{pS} {n1r[x, y]}'
    print(pS)


n2 = s.getPostPen()
print(n2.shape)
n2r = np.round(n2, 2)

for y in range(n2.shape[1]):
    pS = F''
    for x in range(n2.shape[0]):
        pS = F'{pS} {n2r[x, y]}'
    print(pS)

s.calcPostPen(100.0, 0, [0, 10], True, True)

s.printPostPen()
n2 = s.getPostPen()
n2r = np.round(n2, 2)

for y in range(n2.shape[1]):
    pS = F''
    for x in range(n2.shape[0]):
        pS = F'{pS} {n2r[x, y]}'
    print(pS)

print(s.interpolateDistanceImpact(10000, impactIndices.impactVelocity))

print('done')