from pythonwrapper import shell
from pythonwrapper import impactIndices, angleIndices, postPenIndices

import numpy as np
from sklearn.linear_model import LinearRegression

ships = [
    (.100, 1000, .3137,   13, 2154, 10., .010, 17, 45, 60.0, 0, "Akizuki"),
    #(.130, 861, .291,   33.5, 1652, 10,   .03, 22, 45, 60.0, 0, "Okhotnik"),
    (.130, 870, .2857, 33.50, 1700, 10., .010, 22, 45, 60.0, 0, "Leningrad"),
    (.128, 830, .3447,    28, 1640, 10., .010, 21, 45, 60.0, 0, "Maass"),
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
    (.356, 792, .31, 635, 1603, 6, 0.033, 59, 45, 60, 0, "New York")
]
referenceDistancePenetration = {
    "Akizuki":                  {5000:  91, 10000:  47},
    "Okhotnik":                 {5000:  86, 10000:  46},
    "Leningrad":                {5000: 119, 10000:  74, 15000:  47},
    "Maass":                    {5000:  80, 10000:  41},
    "New Orleans":              {5000: 321, 10000: 221, 15000: 154},
    "Budyonny":                 {5000: 215, 10000: 138, 15000:  91},
    "Moskva":                   {5000: 496, 10000: 392, 15000: 314},
    "Nurnberg":                 {5000: 148, 10000:  86, 15000:  52},
    "Graf Spee":                {5000: 405, 10000: 301, 15000: 227},
    "Edinburgh":                {5000: 193, 10000: 118, 15000:  73},
    "Duca d'Aosta":             {5000: 197, 10000: 121, 15000:  76},
    "Arizona":                  {5000: 575, 10000: 467, 15000: 381},
    "North Carolina":           {5000: 662, 10000: 563, 15000: 479},
    "Scharnhorst":              {5000: 455, 10000: 363, 15000: 289},
    "Grosser Kurfurst 420":     {5000: 722, 10000: 624, 15000: 539},
    "Hood":                     {5000: 491, 10000: 406, 15000: 337},
    "King George V":            {5000: 541, 10000: 385, 15000: 326},
    "Thunderer":                {5000: 742, 10000: 646, 15000: 564},
    "Dunkerque":                {5000: 594, 10000: 495, 15000: 414},
    "Oktyabrskaya Revolutsiya": {5000: 458, 10000: 338, 15000: 251},
    "Giulio Cesare":            {5000: 439, 10000: 416, 15000: 320},
    "New York":                 {5000: 337, 15000: 224, 18000: 204},
    #{5000: , 10000:, 15000: },
}

angles = []
velocities = []
diameter = []
mass = []
krupp = []
referenceData = []

def normalization(angle, normalization):
    return max(angle - normalization, 0)


for i, ship in enumerate(ships):
    print(ship[-1])
    s = shell(*ship)
    s.setDtMin(.01)
    s.calcImpactForwardEuler()
    normal = ship[5]
    shipRef = referenceDistancePenetration[ship[-1]]
    for dist, penetration in shipRef.items():
        referenceData.append(penetration)

        adjAngle = normalization(s.interpolateDistanceImpact(dist, int(impactIndices.impactAngleHorizontalDegrees)), normal)
        angles.append(adjAngle)

        impactVelocity = s.interpolateDistanceImpact(dist, int(impactIndices.impactVelocity))
        velocities.append(impactVelocity)
        
        print(impactVelocity, adjAngle)

        diameter.append(ship[0])
        mass.append(ship[3])
        krupp.append(ship[4] / 2400)
    
referenceData = np.array(referenceData)
angles = np.array(angles)
velocities = np.array(velocities)
diameter = np.array(diameter)
mass = np.array(mass)
krupp = np.array(krupp)


def regressV():
    Y = referenceData / np.cos(np.radians(angles)) / krupp / np.power(mass, 0.5506) / np.power(diameter, -0.6521)

    lY = np.log(Y)
    lV = np.log(velocities)

    Xstacked = np.vstack((lV))
    X = Xstacked

    reg = LinearRegression().fit(X, lY)
    print(F'R²: {reg.score(X, lY)}')

    regCoeffs = reg.coef_
    regIntercept = reg.intercept_
    regInterceptE = np.exp(regIntercept)
    print(F'Penetration = {regInterceptE} * V^{regCoeffs[0]} * D^-0.6521 * M^0.5506')

    return regInterceptE * np.cos(np.radians(angles)) * (
                np.power(velocities, regCoeffs[0]) * 
                np.power(diameter, -0.6521) *
                np.power(mass, 0.5506) * 
                np.power(krupp, 1)
            )

def regressVDM():
    Y = referenceData / np.cos(np.radians(angles)) / krupp 

    lY = np.log(Y)
    lV = np.log(velocities)
    lD = np.log(diameter)
    lM = np.log(mass)
    lK = np.log(krupp)

    Xstacked = np.vstack((lV, lD, lM))
    X = Xstacked.T

    reg = LinearRegression().fit(X, lY)
    print(F'R²: {reg.score(X, lY)}')

    regCoeffs = reg.coef_
    regIntercept = reg.intercept_
    regInterceptE = np.exp(regIntercept)
    print(F'Penetration = {regInterceptE} * V^{regCoeffs[0]} * D^{regCoeffs[1]} * M^{regCoeffs[2]}')

    return regInterceptE * np.cos(np.radians(angles)) * (
                np.power(velocities, regCoeffs[0]) * 
                np.power(diameter, regCoeffs[1]) * 
                np.power(mass, regCoeffs[2]) *
                np.power(krupp, 1)
            )

def regressVDMK():
    Y = referenceData / np.cos(np.radians(angles))

    lY = np.log(Y)
    lV = np.log(velocities)
    lD = np.log(diameter)
    lM = np.log(mass)
    lK = np.log(krupp)

    Xstacked = np.vstack((lV, lD, lM, lK))
    X = Xstacked.T

    reg = LinearRegression().fit(X, lY)
    print(F'R²: {reg.score(X, lY)}')

    regCoeffs = reg.coef_
    regIntercept = reg.intercept_
    regInterceptE = np.exp(regIntercept)
    print(F'Penetration = {regInterceptE} * V^{regCoeffs[0]} * D^{regCoeffs[1]} * M^{regCoeffs[2]} * (K/2400)^{regCoeffs[3]}')

    return regInterceptE * np.cos(np.radians(angles)) * (
                np.power(velocities, regCoeffs[0]) * 
                np.power(diameter, regCoeffs[1]) * 
                np.power(mass, regCoeffs[2]) *
                np.power(krupp, regCoeffs[3])
            )

predictions = regressV()

current = 0
for ship in ships:
    row = F'{ship[-1]}: '
    for i in range(len(referenceDistancePenetration[ship[-1]])):
        pCurr = predictions[current]
        row = F'{row} {pCurr} diff: {pCurr - referenceData[current]}; '
        current += 1
    print(row)