from pythonwrapper import shell
from pythonwrapper import impactIndices, angleIndices, postPenIndices

import numpy as np
from sklearn.linear_model import LinearRegression

ships = [
    (.130, 870, .2857, 33.50, 1700, 10., .010, 22, 45, 60.0, 0, "Leningrad"),
    (.203, 853, .3210, 118.0, 2846, 7.0, .033, 34, 60, 67.5, 0, "New Orleans"),
    (.152, 950, .3210, 55.00, 2216, 8.5, .025, 25, 45, 60.0, 0, "Budyonny"),
    (.220, 985, .2549, 176.0, 2590, 7.0, .033, 37, 45, 60.0, 0, "Moskva"),
    (.150, 960, .3307, 45.50, 1862, 8.5, .025, 25, 45, 60.0, 0, "Nurnberg"),
    (.283, 910, .3333, 300.0, 2282, 6.0, .010, 47, 45, 60.0, 0, "Graf Spee"),
    (.152, 841, .3297, 50.8, 2609, 8.5, 0.005, 12, 60, 75.0, 0, "Edinburgh"),
    (.152, 1000, .3256, 50, 2142, 8.5, 0.025, 25, 45, 60.0, 0, "Duca d'Aosta"),
    (.356, 792, .332, 680.4, 2604, 6, 0.033, 59, 45, 60.0, 0, "Arizona"),
    (.406, 701, .352, 1225, 2598, 6, .033, 68, 45, 60.0, 0, "North Carolina"),
    (.283, 890, .2827, 330, 2312, 6, 0.01, 47, 45, 60.0, 0, "Scharnhorst"),
    (.42, 800, .2994, 1220, 2415, 6, .033, 70, 45, 60.0, 0, "Grosser Kurfurst 420"),
    (.381, 731.5, .3379, 879, 2190, 6, 0.033, 64, 45, 60.0, 0, "Hood"),
    #(.356, 757, .3142, 721, 2295, 6, .015, 59, 45, 60.0, 0, "King George V"),
    (.457, 762, .2897, 1506, 2485, 6, .033, 76, 45, 60.0, 0, "Thunderer"),
    (.33, 870, .2801, 560, 2428, 6, .033, 55, 45, 60.0, 0, "Dunkerque"),
    (.305, 762, .4595, 470.9, 2627, 6, 0.01, 51, 45, 60.0, 0, "Oktyabrskaya Revolutsiya"),
    #(.32, 830, .4098, 525, 2600, 6, 0.033, 53, 45, 60.0, 0, "Giulio Cesare"),
]

referenceData = np.array([
    119, 74, 47,   #Leningrad
    321, 221, 154, #New Orleans
    215, 138, 91,  #Budyonny
    496, 392, 314, #Moskva
    148, 86, 52,   #Nurnberg
    405, 301, 227, #Graf Spee
    193, 118, 73,  #Edinburgh
    197, 121, 76,  #Duca d'Aosta
    575, 467, 381, #Arizona
    662, 563, 479, #North Carolina
    455, 363, 289, #Scharnhorst
    722, 624, 539, #Grosser Kurfurst 420
    491, 406, 337, #Hood
    #541, 385, 326, #King George V
    742, 646, 564, #Conqueror - 457
    594, 495, 414, #Dunkerque
    458, 338, 251, #Oktyabrskaya Revolutsiya
    #439, 416, 320  #Guilio Cesare
])

angles = np.zeros(len(ships) * 3)
velocities = np.zeros(len(ships) * 3)
diameter = np.zeros(len(ships) * 3)
mass = np.zeros(len(ships) * 3)
kruppScaling = np.zeros(len(ships) * 3)

def normalization(angle, normalization):
    return max(angle - normalization, 0)


for i, ship in enumerate(ships):
    print(ship[-1])
    s = shell(*ship)
    s.calcImpactForwardEuler()
    #print(s.interpolateDistanceImpact(5000, int(impactIndices.impactVelocity)))
    normal = ship[5]

    angles[3*i  ] = normalization(s.interpolateDistanceImpact( 5000, int(impactIndices.impactAngleHorizontalDegrees)), normal)
    angles[3*i+1] = normalization(s.interpolateDistanceImpact(10000, int(impactIndices.impactAngleHorizontalDegrees)), normal)
    angles[3*i+2] = normalization(s.interpolateDistanceImpact(15000, int(impactIndices.impactAngleHorizontalDegrees)), normal)

    velocities[3*i  ] = s.interpolateDistanceImpact( 5000, int(impactIndices.impactVelocity))
    velocities[3*i+1] = s.interpolateDistanceImpact(10000, int(impactIndices.impactVelocity))
    velocities[3*i+2] = s.interpolateDistanceImpact(15000, int(impactIndices.impactVelocity))

    diameter[3*i] = ship[0]
    diameter[3*i+1] = ship[0]
    diameter[3*i+2] = ship[0]

    mass[3*i] = ship[3]
    mass[3*i+1] = ship[3]
    mass[3*i+2] = ship[3]

    kruppScaling[3*i] = ship[4] / 2400
    kruppScaling[3*i+1] = ship[4] / 2400
    kruppScaling[3*i+2] = ship[4] / 2400

#print(angles, velocities, diameter, mass, kruppScaling)

Y = referenceData / np.cos(np.radians(angles)) / kruppScaling

lY = np.log(Y)
lV = np.log(velocities)
lD = np.log(diameter)
lM = np.log(mass)

X = np.vstack((lV, lD, lM))
print(lV, lD, lM) 
print(X)
print(X.T)
print(lY)


reg = LinearRegression().fit(X.T, lY)
print(reg.score(X.T, lY))

regCoeffs = reg.coef_
regIntercept = reg.intercept_
print(regCoeffs)
print(regIntercept)

print(np.exp(regIntercept))
