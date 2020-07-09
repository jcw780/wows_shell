#ifndef _SHELL_WOWS_HPP_
#define _SHELL_WOWS_HPP_

#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>

#include "controlEnums.hpp"

namespace shell{
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
   public:
    void preProcess() {
        k = 0.5 * cD * pow((caliber / 2), 2) * M_PI / mass;
        // condensed drag coefficient
        cw_2 = 100 + 1000 / 3 * caliber;
        // quadratic drag coefficient
        pPPC = 0.5561613 * krupp / 2400 * pow(mass, 0.5506) /
               pow((caliber * 1000), 0.6521);
        // condensed penetration coefficient
        normalizationR = normalization / 180 * M_PI;
        // normalization (radians)
        ricochet0R = ricochet0 / 180 * M_PI;
        ricochet1R = ricochet1 / 180 * M_PI;
    }

    unsigned int impactSize,
        postPenSize;  // number of distances in: standard, postPen
    unsigned int impactSizeAligned, postPenSizeAligned;
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
    double &get_impact(unsigned int row, unsigned int impact) {
        return impactData[row + impact * impactSizeAligned];
    }

    double *get_impactPtr(unsigned int row, unsigned int impact) {
        return impactData.data() + row + impact * impactSizeAligned;
    }

    double &get_angle(unsigned int row, unsigned int impact) {
        return angleData[row + impact * impactSizeAligned];
    }

    double *get_anglePtr(unsigned int row, unsigned int impact) {
        return angleData.data() + row + impact * impactSizeAligned;
    }

    double &get_postPen(unsigned int row, unsigned int angle,
                        unsigned int impact) {
        return postPenData[row + angle * postPenSize + impact * impactSize];
    }

    double *get_postPenPtr(unsigned int row, unsigned int angle,
                           unsigned int impact) {
        return postPenData.data() + row + angle * postPenSize +
               impact * impactSize;
    }

    // Linear interpolate by distance - WIP

    double interpolateDistanceImpact(double distance, unsigned int impact) {
        auto iter_max = std::lower_bound(
            get_impactPtr(0, impact::impactDataIndex::distance),
            get_impactPtr(0, impact::impactDataIndex::distance + 1), distance);
        double maxDistance = *iter_max;
        unsigned int maxIndex =
            iter_max - get_impactPtr(0, impact::impactDataIndex::distance);
        double maxTarget = get_impact(maxIndex, impact);

        auto iter_min = iter_max - 1;
        double minDistance = *iter_min;
        unsigned int minIndex =
            iter_min - get_impactPtr(0, impact::impactDataIndex::distance);
        double minTarget = get_impact(minIndex, impact);

        /*std::cout << minIndex << " " << minDistance << " " << minTarget << " "
                  << maxIndex << " " << maxDistance << " " << maxTarget << "\n";
                  */
        double slope = ((maxTarget - minTarget) / (maxDistance - minDistance));
        return slope * (distance - minDistance) + minTarget;
    }

    // internal computed data - fixed
    double get_v0() { return v0; }
    double get_k() { return k; }
    double get_cw_2() { return cw_2; }
    double get_pPPC() { return pPPC; }
    double get_normalizationR() { return normalizationR; }

    void printAngleData() {
        for (unsigned int i = 0; i < impactSize; i++) {
            for (unsigned int j = 0; j < angle::maxColumns; j++) {
                std::cout << std::fixed << std::setprecision(4)
                          << get_angle(i, j) << " ";
            }
            std::cout << "\n";
        }
        std::cout << "Completed Angle Data\n";
    }

    void printPostPenData() {
        for (unsigned int i = 0; i < postPenSize; i++) {
            for (unsigned int j = 0; j < post::maxColumns; j++) {
                std::cout << std::fixed << std::setprecision(4)
                          << get_postPen(i, j, 0) << " ";
            }
            std::cout << "\n";
        }
        std::cout << "Completed Post-Penetration\n";
    }

    void printImpactData() {
        for (unsigned int i = 0; i < impactSize; i++) {
            for (unsigned int j = 0; j < impact::maxColumns; j++) {
                std::cout << std::fixed << std::setprecision(4)
                          << get_impact(i, j) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "Completed Standard Data" << std::endl;
    }
    void printTrajectory(unsigned int target) {
        if (target >= impactSize) {
            std::cout << "Target Not Within Range of: " << impactSize
                      << std::endl;
        } else {
            printf("Index:[%d] X Y\n", target);
            for (std::vector<double>::size_type i = 0;
                 i < trajectories[target * 2].size(); i++) {
                std::cout << trajectories[target * 2][i] << " "
                          << trajectories[target * 2 + 1][i] << std::endl;
            }
        }
    }
};
}
#endif