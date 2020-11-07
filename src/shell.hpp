#pragma once

#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

#include "controlEnums.hpp"
#include "utility.hpp"

namespace wows_shell {
double combinedAirDrag(double cD, double caliber, double mass){
    return 0.5 * cD * pow((caliber / 2), 2) * M_PI / mass;
}

double combinedPenetration(double krupp, double mass, double caliber){
    return 0.00046905491615181766 * krupp / 2400 * pow(mass, 0.5506) *
               pow(caliber, -0.6521);
}

class shell {
   public:                 // Description                units
    double v0;             // muzzle velocity            m/s
    double caliber;        // shell caliber              m
    double krupp;          // shell krupp                [ndim]
    double mass;           // shell mass                 kg
    double cD;             // coefficient of drag        [ndim]
    double normalization;  // shell normalization        degrees
    double threshold;      // fusing threshold           mm
    double fuseTime;       // fuse time                  s
    double ricochet0;      // ricochet angle 0           degrees
    double ricochet1;      // ricochet angle 1           degrees
    double nonAP;          // penetration of non ap ammo mm
    bool enableNonAP;
    std::string name;

    double k, cw_2, pPPC, normalizationR, ricochet0R, ricochet1R;

    // Condenses initial values into values used by calculations
    //[Reduces repeated computations]
    void preProcess() {
        k = combinedAirDrag(cD, caliber, mass);
        //k = 0.5 * cD * pow((caliber / 2), 2) * M_PI / mass;
        // condensed drag coefficient
        // cw_2 = 100 + 1000 / 3 * caliber;
        cw_2 = 0;
        // linear drag coefficient
        // pPPC = 0.5561613 * krupp / 2400 * pow(mass, 0.5506) /
        //       pow((caliber * 1000), 0.6521);
        // pPPC = 0.081525 * krupp / 2400 * pow(mass, 0.5506) /
        //       pow((caliber * 1000), 0.6521);
        //pPPC = 0.00046905491615181766 * krupp / 2400 * pow(mass, 0.5506) *
        //       pow(caliber, -0.6521);
        pPPC = combinedAirDrag(krupp, caliber, mass);

        // condensed penetration coefficient
        normalizationR = normalization / 180 * M_PI;
        // normalization (radians)
        ricochet0R = ricochet0 / 180 * M_PI;
        ricochet1R = ricochet1 / 180 * M_PI;
    }

    // Not 100% necessary - sizes adjusted to fulfill alignment
    bool completedImpact = false, completedAngles = false,
         completedPostPen = false;

    /*trajectories output
    [0           ]trajx 0        [1           ]trajy 1
    ...
    [size * 2 - 2]trajx size - 1 [size * 2 - 1]trajy size - 1
    */
    std::vector<std::vector<double>> trajectories;

    // Refer to stdDataIndex enums defined above
    std::size_t impactSize = 0, impactSizeAligned;
    std::vector<double> impactData;

    /* Angles data
     * [0:1)-ra0 max lateral angle
     * [1:2)-ra0 max lateral angle
     * [2:3)-ra1 max lateral angle
     * [3:4)-ra1 max lateral angle
     * [4:5)-penetration max lateral angle
     * [5:6)-penetration max lateral angle
     */
    std::vector<double> angleData;

    /* WARNING: LOCATION OF LATERAL ANGLE IN VECTOR CANNOT BE CHANGED OR ELSE
     * SIMD ALIGNMENT MAY NOT BE GUARANTEED [0:1) Lateral Angle [1:2) Distance
     * [2:3) X [3:4) Y [4:5) Z             [5:6) XWF See enums defined above
     */
    std::size_t postPenSize = 0, postPenSizeAligned;
    std::vector<double> postPenData;

    shell() = default;

    shell(const double caliber, const double v0, const double cD,
          const double mass, const double krupp, const double normalization,
          const double fuseTime, const double threshold, const double ricochet0,
          const double ricochet1, const double nonAP, const std::string &name) {
        setValues(caliber, v0, cD, mass, krupp, normalization, fuseTime,
                  threshold, ricochet0, ricochet1, nonAP, name);
    }

    void setValues(const double caliber, const double v0, const double cD,
                   const double mass, const double krupp,
                   const double normalization, const double fuseTime,
                   const double threshold, const double ricochet0,
                   const double ricochet1, const double nonAP,
                   const std::string &name) {
        //                                                    Impact   PostPen
        this->fuseTime = fuseTime;  // Shell fusetime        | No     | Yes
        this->v0 = v0;              // Shell muzzle velocity | Yes    | No
        this->caliber = caliber;    // Shell caliber         | Yes    | Yes
        this->krupp = krupp;        // Shell krupp           | Yes    | Yes
        this->mass = mass;          // Shell mass            | Yes    | Yes
        this->normalization =
            normalization;  // Shell normalization   | Yes    | Yes
        this->cD = cD;      // Shell air drag coefficient | Yes    | Yes
        this->name = name;  // Shell name  | No     | No
        this->ricochet0 = ricochet0;  // Ricochet Angle 0 | No     | Yes
        this->ricochet1 = ricochet1;  // Ricochet Angle 1 | No     | Yes
                                      // Shell fusing threshold | No     | Yes
        this->threshold = threshold;
        this->nonAP = nonAP;
        if (this->nonAP > 0) {
            this->enableNonAP = true;
        } else {
            this->enableNonAP = false;
        }
        preProcess();
    }

    // Getter Functions
    double &get_impact(const std::size_t row, impact::impactIndices impact) {
        return get_impact(row, toUnderlying(impact));
    }
    double &get_impact(const std::size_t row, const std::size_t impact) {
        return impactData[row + impact * impactSizeAligned];
    }

    double *get_impactPtr(const std::size_t row, impact::impactIndices impact) {
        return get_impactPtr(row, toUnderlying(impact));
    }
    double *get_impactPtr(const std::size_t row, const std::size_t impact) {
        return impactData.data() + row + impact * impactSizeAligned;
    }

    double &get_angle(const std::size_t row, angle::angleIndices data) {
        return get_angle(row, toUnderlying(data));
    }
    double &get_angle(const std::size_t row, const std::size_t impact) {
        return angleData[row + impact * impactSizeAligned];
    }

    double *get_anglePtr(const std::size_t row, angle::angleIndices data) {
        return get_anglePtr(row, toUnderlying(data));
    }
    double *get_anglePtr(const std::size_t row, const std::size_t impact) {
        return angleData.data() + row + impact * impactSizeAligned;
    }

    double &get_postPen(const std::size_t row, post::postPenIndices data,
                        const std::size_t angle) {
        return get_postPen(row, toUnderlying(data), angle);
    }
    double &get_postPen(const std::size_t row, const std::size_t data,
                        const std::size_t angle) {
        return postPenData[row + data * postPenSize + angle * impactSize];
    }

    double *get_postPenPtr(const std::size_t row, post::postPenIndices data,
                           const std::size_t angle) {
        return get_postPenPtr(row, toUnderlying(data), angle);
    }
    double *get_postPenPtr(const std::size_t row, const std::size_t angle,
                           const std::size_t impact) {
        return postPenData.data() + row + angle * postPenSize +
               impact * impactSize;
    }

    std::size_t maxDist() {
        std::size_t errorCode = std::numeric_limits<std::size_t>::max();
        if (impactSize == 0) return errorCode;
        if (impactSize == 1) return 0;
        if (impactSize == 2)
            return get_impact(1, impact::impactIndices::distance) >
                           get_impact(0, impact::impactIndices::distance)
                       ? 1
                       : 0;

        std::size_t start = 0, end = impactSize - 1;
        if (get_impact(end, impact::impactIndices::distance) >=
            get_impact(end - 1, impact::impactIndices::distance))
            return end;

        using T = typename std::decay<decltype(*impactData.begin())>::type;
        while (start <= end) {
            std::size_t mid = start + (end - start) / 2;
            T midValue = get_impact(mid, impact::impactIndices::distance),
              leftValue = get_impact(mid - 1, impact::impactIndices::distance),
              rightValue = get_impact(mid + 1, impact::impactIndices::distance);
            // std::cout << mid << " " << leftValue << " " << midValue << " "
            //          << rightValue << "\n";
            if (leftValue <= midValue && midValue >= rightValue) {
                return mid;
            }
            if (leftValue < midValue && midValue < rightValue) {
                start = ++mid;
            } else if (leftValue > midValue && midValue > rightValue) {
                end = --mid;
            }
        }
        return errorCode;
    }

    double interpolateDistanceImpact(double distance,
                                     impact::impactIndices data) {
        return interpolateDistanceImpact(distance, toUnderlying(data));
    }
    double interpolateDistanceImpact(double distance, uint32_t impact) {
        std::size_t maxIndex = maxDist(),
                    maxErrorCode = std::numeric_limits<std::size_t>::max();
        double errorCode = std::numeric_limits<double>::max();
        if (maxIndex == maxErrorCode) return errorCode;
        if (distance < get_impact(0, impact::impactIndices::distance))
            return errorCode;
        if (distance > get_impact(maxIndex, impact::impactIndices::distance))
            return errorCode;

        // Only get lower set
        auto iter_max = std::lower_bound(
            get_impactPtr(0, impact::impactIndices::distance),
            get_impactPtr(maxIndex, impact::impactIndices::distance), distance);
        double upperDistance = *iter_max;
        uint32_t upperIndex =
            iter_max - get_impactPtr(0, impact::impactIndices::distance);
        double upperTarget = get_impact(upperIndex, impact);

        if (upperIndex == 0) return upperIndex;
        // Only activates if distance = min and prevents segfault

        auto iter_min = iter_max - 1;
        double lowerDistance = *iter_min;
        uint32_t lowerIndex =
            iter_min - get_impactPtr(0, impact::impactIndices::distance);
        double lowerTarget = get_impact(lowerIndex, impact);

        /*std::cout << minIndex << " " << minDistance << " " << minTarget << " "
                  << maxIndex << " " << maxDistance << " " << maxTarget << "\n";
                  */
        double slope =
            ((upperTarget - lowerTarget) / (upperDistance - lowerDistance));
        return slope * (distance - lowerDistance) + lowerTarget;
    }

    // internal computed data - fixed
    const double &get_v0() { return v0; }
    const double &get_k() { return k; }
    const double &get_cw_2() { return cw_2; }
    const double &get_pPPC() { return pPPC; }
    const double &get_normalizationR() { return normalizationR; }

    void printAngleData() {
        for (std::size_t i = 0; i < impactSize; i++) {
            for (std::size_t j = 0; j < angle::maxColumns; j++) {
                std::cout << std::fixed << std::setprecision(4)
                          << get_angle(i, j) << " ";
            }
            std::cout << "\n";
        }
        std::cout << "Completed Angle Data\n";
    }

    void printPostPenData() {
        for (std::size_t i = 0; i < postPenSize; i++) {
            for (std::size_t j = 0; j < post::maxColumns; j++) {
                std::cout << std::fixed << std::setprecision(4)
                          << get_postPen(i, j, 0) << " ";
            }
            std::cout << "\n";
        }
        std::cout << "Completed Post-Penetration\n";
    }

    void printImpactData() {
        for (std::size_t i = 0; i < impactSize; i++) {
            for (std::size_t j = 0; j < impact::maxColumns; j++) {
                std::cout << std::fixed << std::setprecision(4)
                          << get_impact(i, j) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "Completed Standard Data" << std::endl;
    }
    void printTrajectory(std::size_t target) {
        if (target >= impactSize) {
            std::cout << "Target Not Within Range of: " << impactSize
                      << std::endl;
        } else {
            std::cout << "Index:[" << target << "] X Y\n";
            for (std::vector<double>::size_type i = 0;
                 i < trajectories[target * 2].size(); i++) {
                std::cout << trajectories[target * 2][i] << " "
                          << trajectories[target * 2 + 1][i] << std::endl;
            }
        }
    }
};

std::string generateHash(const double k, const double p, const double v0, const double normalization,
          const double fuseTime, const double threshold, const double ricochet0,
          const double ricochet1, const double nonAP){
    std::vector<char> hashString(8 * sizeof(double) + sizeof(bool), 0);
    memcpy(&hashString[0*sizeof(double)], &k, sizeof(double));
    memcpy(&hashString[1*sizeof(double)], &p, sizeof(double));
    memcpy(&hashString[2*sizeof(double)], &v0, sizeof(double));
    memcpy(&hashString[3*sizeof(double)], &normalization, sizeof(double));
    memcpy(&hashString[4*sizeof(double)], &threshold, sizeof(double));
    memcpy(&hashString[5*sizeof(double)], &fuseTime, sizeof(double));
    memcpy(&hashString[6*sizeof(double)], &ricochet0, sizeof(double));
    memcpy(&hashString[7*sizeof(double)], &ricochet1, sizeof(double));
    memcpy(&hashString[8*sizeof(double)], &nonAP, sizeof(bool));
    
    return utility::base64_encode(hashString);
}

std::string generateHash(const double caliber, const double v0, const double cD,
          const double mass, const double krupp, const double normalization,
          const double fuseTime, const double threshold, const double ricochet0,
          const double ricochet1, const double nonAP){
    double k = combinedAirDrag(cD, caliber, mass);
    double p = combinedPenetration(krupp, mass, caliber);
    return generateHash(k, p, v0, normalization, threshold, fuseTime, ricochet0, ricochet1, nonAP);
}

std::string generateHash(const shell& s){
    return generateHash(
        s.caliber, s.v0, s.cD, 
        s.mass, s.krupp, s.normalization, 
        s.fuseTime, s.threshold, s.ricochet0, 
        s.ricochet1, s.nonAP);
}

}  // namespace wows_shell