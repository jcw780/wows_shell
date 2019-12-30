from pythonwrapper import shell
import numpy as np
s = shell(780, .460, 2574, 1460, 6, .292, "Yamato", 76.0, .033 )
s.calcImpact()
print("Standard Done")
s.calcPostPen(400.0, [0, 10, 20])
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

s.calcPostPen(100.0, [0, 10])

s.printPostPen()
n2 = s.getPostPen()
n2r = np.round(n2, 2)

for y in range(n2.shape[1]):
    pS = F''
    for x in range(n2.shape[0]):
        pS = F'{pS} {n2r[x, y]}'
    print(pS)

print('done')