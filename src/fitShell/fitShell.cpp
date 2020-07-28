#include "../shellCPP.hpp"

namespace shell {
enum sampleIndices { launchA, distance };
namespace fitPenetration {
enum index { distance, penetration };
}
inline int signum(double x) { return ((0.0) < x) - (x < (0.0)); }

template <unsigned int Numerical>
void fitDrag(shell &toFit, std::vector<double> &sampleData,
             unsigned int length) {
    shellCalc calculator;
    calculator.set_dt_min(.0001);
    toFit.impactSize = length;
    toFit.impactSizeAligned = calculator.calculateAlignmentSize(length);
    toFit.impactData.resize(toFit.impactSizeAligned * impact::maxColumnsFit);
    std::copy_n(&sampleData[sampleIndices::launchA], length,
                toFit.get_impactPtr(0, impact::impactIndices::launchAngle));

    double currentCD = .35;
    double learningRate = .000004;

    int sign = 1;

    double gradient = calculateGradient<Numerical>(toFit, sampleData, length,
                                                   calculator, currentCD);
    sign = signum(gradient);
    std::cout << std::setprecision(10) << 0 << " cD: " << currentCD
              << " gradient: " << gradient << "\n";
    currentCD -= learningRate * gradient;

    for (int i = 1; i < 30; i++) {
        double gradient = calculateGradient<Numerical>(
            toFit, sampleData, length, calculator, currentCD);
        std::cout << std::setprecision(10) << i << " cD: " << currentCD
                  << " gradient: " << gradient << "\n";
        int gradientSign = signum(gradient);
        if (gradientSign != sign) {
            learningRate /= 2;
            sign = gradientSign;
        }

        currentCD -= learningRate * gradient;
    }
    toFit.cD = currentCD;
    toFit.preProcess();
}
template <unsigned int Numerical>
double calculateGradient(shell &toFit, std::vector<double> &sampleData,
                         unsigned int length, shellCalc &calculator,
                         double currentCD) {
    std::array<double, 2> errors;
    static constexpr double deltaCD = .001;
    toFit.cD = (currentCD - deltaCD);
    toFit.preProcess();
    calculator.calculateFit<Numerical>(toFit);
    errors[0] = calculateErrors(toFit, sampleData, length);
    toFit.cD = (currentCD + deltaCD);
    toFit.preProcess();
    calculator.calculateFit<Numerical>(toFit);
    errors[1] = calculateErrors(toFit, sampleData, length);
    std::cout << "stddev: " << errors[0] << " " << errors[1] << "\n";
    return (errors[1] - errors[0]) / (2 * deltaCD);
}

double calculateErrors(shell &toFit, std::vector<double> &sampleData,
                       unsigned int length) {
    double errors = 0;
    for (unsigned int i = 0; i < length; i++) {
        errors += pow(toFit.get_impact(i, impact::impactIndices::distance) -
                          sampleData[i + length * sampleIndices::distance],
                      2);
    }
    return sqrt(errors / length);
}

double normalizationCos(double angleRadians, double normalizationDegrees) {
    double normalizationR = normalizationDegrees * M_PI / 180;
    double returnV = (fabs(angleRadians) > normalizationR) *
                     (fabs(angleRadians) - normalizationR);
    // std::cout << "NormalV: " << returnV << " " << normalizationR << " ";
    return cos(returnV);
}

void gradientKruppNormal(const std::vector<double> &velocityAngleData,
                         const std::vector<double> &sampleData,
                         const unsigned int length, const double &krupp,
                         const double &pPPCmK, const double &normal,
                         std::array<double, 2> &gradients) {
    std::array<double, 4> errors{0, 0, 0, 0};
    static constexpr double dKrupp = 1;
    std::array<double, 2> krupps = {krupp - dKrupp, krupp + dKrupp};
    static constexpr double dNormal = .5;
    std::array<double, 2> normals = {normal - dNormal, normal + dNormal};
    for (unsigned int i = 0; i < length; i++) {
        const double fallAngle = velocityAngleData[i + length];
        const double referencePenetration =
            sampleData[fitPenetration::index::penetration * length + i];
        const double eAR = normalizationCos(fallAngle, normal);
        const double velocityPPCMK = velocityAngleData[i] * pPPCmK;
        const double penetrationCoeff = velocityPPCMK * krupp;

        for (int j = 0; j < 2; j++) {
            errors[j] += pow(
                (referencePenetration - (velocityPPCMK * krupps[j] * eAR)), 2);
            errors[j + 2] += pow(
                (referencePenetration -
                 (penetrationCoeff * normalizationCos(fallAngle, normals[j]))),
                2);
        }
        // std::cout << errors[0] << " " << errors[1] << " " << errors[2] << " "
        //          << errors[3] << "\n";
    }
    std::cout << "stddev: ";
    for (int i = 0; i < 4; i++) {
        errors[i] = sqrt(errors[i] / length);
        std::cout << errors[i] << " ";
    }
    std::cout << "nl\n";
    gradients[0] = (errors[1] - errors[0]) / (2 * dKrupp);
    gradients[1] = (errors[3] - errors[2]) / (2 * dNormal);
    std::cout << "Krupp: " << krupp << " Gradient: " << gradients[0]
              << " Normal: " << normal << " Gradient: " << gradients[1] << "\n";
}

void fitKruppNormal(shell &toFit, std::vector<double> &sampleData,
                    unsigned int length, double maxAngle) {
    shellCalc calculator;
    calculator.set_max(maxAngle);
    calculator.set_precision(.1);
    calculator.calculateImpact<false, numerical::rungeKutta2, false>(toFit);
    std::vector<double> velocityAngleData(length * 2);
    for (unsigned int i = 0; i < length; i++) {
        double tgtDist =
            sampleData[fitPenetration::index::distance * length + i];
        velocityAngleData[i] =
            pow(toFit.interpolateDistanceImpact(
                    tgtDist, impact::impactIndices::impactVelocity),
                1.1);
        velocityAngleData[i + length] = toFit.interpolateDistanceImpact(
            tgtDist, impact::impactIndices::impactAngleHorizontalRadians);
        std::cout << tgtDist << " " << velocityAngleData[i] << " "
                  << velocityAngleData[i + length] << "\n";
    }

    std::array<double, 2> learningRates{2000, .1};
    double krupp = 2400;
    double normalization = 10;

    const double pPPCmK = 0.5561613 / 2400 * pow(toFit.mass, 0.55) /
                          pow((toFit.caliber * 1000), 0.65);
    std::cout << pPPCmK << "pPPCmK \n";

    std::array<int, 2> signs;
    std::array<double, 2> gradients;
    {
        std::cout << 1 << "\n";
        gradientKruppNormal(velocityAngleData, sampleData, length, krupp,
                            pPPCmK, normalization, gradients);
        for (int j = 0; j < 2; j++) {
            signs[j] = calculator.signum(gradients[j]);
        }
        krupp -= learningRates[0] * gradients[0];
        normalization -= learningRates[1] * gradients[1];
        std::cout << "learningRates: " << learningRates[0] << " "
                  << learningRates[1] << "\n";
        std::cout << "LRgradients: " << -1 * learningRates[0] * gradients[0]
                  << " " << -1 * learningRates[1] * gradients[1] << "\n";
    }

    for (int i = 1; i < 200; i++) {
        std::cout << i + 1 << "\n";
        gradientKruppNormal(velocityAngleData, sampleData, length, krupp,
                            pPPCmK, normalization, gradients);
        for (int j = 0; j < 2; j++) {
            int signsT = calculator.signum(gradients[j]);
            if (signs[j] != signsT) {
                learningRates[j] /= 2;
                signs[j] = signsT;
            }
        }
        krupp -= learningRates[0] * gradients[0];
        normalization -= learningRates[1] * gradients[1];
        std::cout << "learningRates: " << learningRates[0] << " "
                  << learningRates[1] << "\n";
        std::cout << "deltas: " << -1 * learningRates[0] * gradients[0] << " "
                  << -1 * learningRates[1] * gradients[1] << "\n";
    }

    std::cout << "final krupp: " << krupp
              << " final normalization: " << normalization << "\n";
    toFit.krupp = krupp;
    toFit.normalization = normalization;
    toFit.preProcess();
}

}  // namespace shell

int main() {
    shell::shell test(.406, 762, 0.2988329237, 1225, 2520, 6, .033, 76, 45, 60,
                      0, "Montana");
    std::vector<double> sample = {10,    15,    20,    25,    30,    35,
                                  40,    45,    16139, 21854, 26518, 30450,
                                  33558, 36119, 37884, 38720};
    shell::fitDrag<shell::numerical::adamsBashforth5>(test, sample, 8);
    std::cout << "cD: " << std::setprecision(10) << test.cD << "\n";

    /*std::vector<double> penetrationData = {
        4572, 9144, 13716, 18288, 22860, 27432, 32004, 36576, 38720,
        747,  664,  585,   509,   441,   380,   329,   280,   241};

    shell::fitKruppNormal(test, penetrationData, 9, 50);*/

    shell::shellCalc scV;
    scV.set_precision(5);
    scV.set_max(45.0);
    scV.calculateImpact<false, shell::numerical::rungeKutta2, false>(test);
    test.printImpactData();
    return 0;
}