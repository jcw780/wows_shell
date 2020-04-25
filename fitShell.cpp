#include "shellCPP.hpp"

namespace shell {
enum sampleDataIndex { launchA, distance };
namespace fitPenetration {
enum fitPenetrationIndex { distance, penetration };
}

template <unsigned int Numerical>
void fitDrag(shell &toFit, std::vector<double> &sampleData,
             unsigned int length) {
    shellCalc calculator;
    calculator.set_dt_min(.0001);
    toFit.impactSize = length;
    toFit.impactSizeAligned = calculator.calculateAlignmentSize(length);
    toFit.impactData.resize(toFit.impactSizeAligned * impact::maxColumnsFit);
    std::copy_n(&sampleData[sampleDataIndex::launchA], length,
                toFit.get_impactPtr(0, impact::impactDataIndex::launchA));

    double currentCD = .35;
    double learningRate = .000004;

    int sign = 1;

    double gradient = calculateGradient<Numerical>(toFit, sampleData, length,
                                                   calculator, currentCD);
    sign = calculator.signum(gradient);
    std::cout << std::setprecision(10) << 0 << " cD: " << currentCD
              << " gradient: " << gradient << "\n";
    currentCD -= learningRate * gradient;

    for (int i = 1; i < 30; i++) {
        double gradient = calculateGradient<Numerical>(
            toFit, sampleData, length, calculator, currentCD);
        std::cout << std::setprecision(10) << i << " cD: " << currentCD
                  << " gradient: " << gradient << "\n";
        int gradientSign = calculator.signum(gradient);
        if (gradientSign != sign) {
            learningRate /= 2;
            sign = gradientSign;
        }

        currentCD -= learningRate * gradient;
    }
    toFit.set_cD(currentCD);
    toFit.preProcess();
}
template <unsigned int Numerical>
double calculateGradient(shell &toFit, std::vector<double> &sampleData,
                         unsigned int length, shellCalc &calculator,
                         double currentCD) {
    std::array<double, 2> errors;
    static constexpr double deltaCD = .001;
    toFit.set_cD(currentCD - deltaCD);
    toFit.preProcess();
    calculator.calculateFit<Numerical>(toFit);
    errors[0] = calculateErrors(toFit, sampleData, length);
    toFit.set_cD(currentCD + deltaCD);
    toFit.preProcess();
    calculator.calculateFit<Numerical>(toFit);
    errors[1] = calculateErrors(toFit, sampleData, length);
    std::cout << "stddev: " << errors[0] << " " << errors[1] << "\n";
    return (errors[1] - errors[0]) / (2 * deltaCD);
}

double calculateErrors(shell &toFit, std::vector<double> &sampleData,
                       unsigned int length) {
    double errors;
    for (int i = 0; i < length; i++) {
        errors += pow(toFit.get_impact(i, impact::impactDataIndex::distance) -
                          sampleData[i + length * sampleDataIndex::distance],
                      2);
    }
    return sqrt(errors / length);
}

void fitPenetration(shell &toFit, std::vector<double> &sampleData,
                    unsigned int length, double maxAngle) {
    double generalPenetrationCoefficient;
    double normalization;
    shellCalc calculator;
    calculator.set_max(maxAngle);
    std::vector<double> velocities(length);
    for (int i = 0; i < length; i++) {
        double distance;
    }
}

} // namespace shell

int main() {
    shell::shell test(.406, 762, .292, 1225, 2520, 6, .033, 76, 45, 60,
                      "Montana");
    std::vector<double> sample = {10,    15,    20,    25,    30,    35,
                                  40,    45,    16139, 21854, 26518, 30450,
                                  33558, 36119, 37884, 38720};
    shell::fitDrag<shell::numerical::rungeKutta2>(test, sample, 8);
    std::cout << std::setprecision(10) << test.cD << "\n";

    shell::shellCalc scV;
    scV.set_precision(5);
    scV.set_max(46.0);
    scV.calculateImpact<false, shell::numerical::rungeKutta2, false>(test);
    test.printImpactData();
    return 0;
}