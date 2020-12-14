#pragma once

#define _USE_MATH_DEFINES
#include <algorithm>
#include <array>
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
struct shellParams {
    double caliber;
    double v0;
    double cD;
    double mass;
    double krupp;
    double normalization;
    double fuseTime;
    double threshold;
    double ricochet0;
    double ricochet1;
    double nonAP;
};

struct dispersionParams {
    double idealRadius;    // maxRadius at ideal distance 30m
    double minRadius;      // minRadius at dist: 0 30m
    double idealDistance;  // idealDistance 30m
    double taperDistance;  // taperDistance m
    double delim;          // vertical seperator ratio [1]
    double zeroRadius;     // vertical ratio at dist: 0 [1]
    double delimRadius;    // vertical ratio at delim [1]
    double maxRadius;      // vertical ratio at maxDistance [1]
    double maxDistance;    // max distance [m]
    double sigma;          // bound of truncated normal distribution
};

double combinedAirDrag(double cD, double caliber, double mass) {
    return 0.5 * cD * pow((caliber / 2), 2) * M_PI / mass;
}

double combinedPenetration(double krupp, double mass, double caliber) {
    return 0.00046905491615181766 * krupp / 2400 * pow(mass, 0.5506) *
           pow(caliber, -0.6521);
}

class shell {
   public:  // Description                units
    // Shell
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

    // Dispersion
    // Horizontal
    double idealRadius;
    double minRadius;
    double idealDistance;
    double taperDistance;
    // Intermediate Horizontal Values
    double taperSlope;
    double horizontalSlope;
    double horizontalIntercept;

    // Vertical
    double delim;
    double zeroRadius;
    double delimRadius;
    double maxRadius;
    double maxDistance;
    // Intermediate Vertical Values
    double zeroDelimSlope;
    double zeroDelimIntercept;
    double delimMaxSlope;
    double delimMaxIntercept;
    double delimDistance;
    bool convex;

    // Distribution
    double sigma;
    // Intermediate Sigma Values
    double standardRatio;
    double halfRatio;

    bool enableNonAP;
    std::string name;

    double k, cw_2, pPPC, normalizationR, ricochet0R, ricochet1R;

    // Condenses initial values into values used by calculations
    //[Reduces repeated computations]
    void preProcess() {
        k = combinedAirDrag(cD, caliber, mass);  // condensed drag coefficient
        cw_2 = 0;                                // linear drag coefficient
        pPPC = combinedPenetration(krupp, mass, caliber);
        // condensed penetration coefficient
        normalizationR = normalization / 180 * M_PI;
        // normalization (radians)
        ricochet0R = ricochet0 / 180 * M_PI;
        ricochet1R = ricochet1 / 180 * M_PI;

        horizontalSlope = (idealRadius - minRadius) / idealDistance;
        horizontalIntercept = 30 * minRadius;
        taperSlope = (horizontalSlope * taperDistance + horizontalIntercept) /
                     taperDistance;

        delimDistance = maxDistance * delim;
        zeroDelimSlope = (delimRadius - zeroRadius) / delimDistance;
        zeroDelimIntercept = zeroRadius;
        delimMaxSlope = (maxRadius - delimRadius) / delimDistance;
        delimMaxIntercept = delimRadius - delimMaxSlope * delimDistance;
        convex = zeroDelimSlope >= delimMaxSlope;
        // This allows for an optimization that removes branches in a hot loop

        const double left = -1 * sigma, right = sigma;
        const double phiL = utility::pdf(left), phiR = utility::pdf(right);
        const double Z = utility::cdf(right) - utility::cdf(left);
        const double DsQ = (phiL - phiR) / Z;
        standardRatio =
            sqrt(1 + ((left * phiL - right * phiR) / Z) - (DsQ * DsQ)) / sigma;
        halfRatio = utility::invCDF(.25 * Z + phiL) * -1;
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

    /* Dispersion data
     *  [0:1)- max horizontal
     *  [1:2)- std horizontal
     *  [2:3)- 50% horizontal
     *  [3:4)- max vertical
     *  [4:5)- std vertical
     *  [5:6)- 50% vertical
     *  [6:7)- max area
     *  [7:8)- std area
     *  [8:9)- 50% area
     */
    std::vector<double> dispersionData;

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
    }  // TODO: Deprecate / Remove because this is very unclean

    shell(const shellParams sp, const std::string &name) {
        setValues(sp, name);
    }

    shell(const shellParams sp, const dispersionParams dp,
          const std::string &name) {
        setValues(sp, dp, name);
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

    void setValues(const shellParams sp, const std::string &name) {
        fuseTime = sp.fuseTime;  // Shell fusetime        | No     | Yes
        v0 = sp.v0;              // Shell muzzle velocity | Yes    | No
        caliber = sp.caliber;    // Shell caliber         | Yes    | Yes
        krupp = sp.krupp;        // Shell krupp           | Yes    | Yes
        mass = sp.mass;          // Shell mass            | Yes    | Yes
        normalization =
            sp.normalization;      // Shell normalization   | Yes    | Yes
        cD = sp.cD;                // Shell air drag coefficient | Yes    | Yes
        this->name = name;         // Shell name  | No     | No
        ricochet0 = sp.ricochet0;  // Ricochet Angle 0 | No     | Yes
        ricochet1 = sp.ricochet1;  // Ricochet Angle 1 | No     | Yes
                                   // Shell fusing threshold | No     | Yes
        threshold = sp.threshold;
        nonAP = sp.nonAP;
        if (nonAP > 0) {
            enableNonAP = true;
        } else {
            enableNonAP = false;
        }
        preProcess();
    }

    void setValues(const shellParams sp, const dispersionParams dp,
                   const std::string &name) {
        fuseTime = sp.fuseTime;  // Shell fusetime        | No     | Yes
        v0 = sp.v0;              // Shell muzzle velocity | Yes    | No
        caliber = sp.caliber;    // Shell caliber         | Yes    | Yes
        krupp = sp.krupp;        // Shell krupp           | Yes    | Yes
        mass = sp.mass;          // Shell mass            | Yes    | Yes
        normalization =
            sp.normalization;      // Shell normalization   | Yes    | Yes
        cD = sp.cD;                // Shell air drag coefficient | Yes    | Yes
        this->name = name;         // Shell name  | No     | No
        ricochet0 = sp.ricochet0;  // Ricochet Angle 0 | No     | Yes
        ricochet1 = sp.ricochet1;  // Ricochet Angle 1 | No     | Yes
                                   // Shell fusing threshold | No     | Yes
        threshold = sp.threshold;
        nonAP = sp.nonAP;
        if (nonAP > 0) {
            enableNonAP = true;
        } else {
            enableNonAP = false;
        }

        idealRadius = dp.idealRadius;
        minRadius = dp.minRadius;
        idealDistance = dp.idealDistance;
        taperDistance = dp.taperDistance;
        delim = dp.delim;
        zeroRadius = dp.zeroRadius;
        delimRadius = dp.delimRadius;
        maxRadius = dp.maxRadius;
        maxDistance = dp.maxDistance;
        sigma = dp.sigma;

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

    double *get_dispersionPtr(const std::size_t row, const std::size_t impact) {
        return dispersionData.data() + row + impact * impactSizeAligned;
    }

    double *get_dispersionPtr(const std::size_t row,
                              dispersion::dispersionIndices data) {
        return get_dispersionPtr(row, toUnderlying(data));
    }

    double &get_dispersion(const std::size_t row, const std::size_t impact) {
        return *get_dispersionPtr(row, impact);
    }

    double &get_dispersion(const std::size_t row,
                           dispersion::dispersionIndices data) {
        return *get_dispersionPtr(row, data);
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

    void printDispersionData() {
        for (std::size_t i = 0; i < impactSize; ++i) {
            std::cout << std::fixed << std::setprecision(4)
                      << get_impact(i, impact::impactIndices::distance) << " ";
            for (std::size_t j = 0; j < dispersion::maxColumns; ++j) {
                std::cout << std::fixed << std::setprecision(4)
                          << get_dispersion(i, j) << " ";
            }
            std::cout << "\n";
        }
        std::cout << "Completed Dispersion Data\n";
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

template <typename T>
std::string generateHash(const T k, const T p, const T v0,
                         const T normalization, const T fuseTime,
                         const T threshold, const T ricochet0,
                         const T ricochet1, const T nonAP) {
    std::array<char, 9 * sizeof(T)> hashString;
    std::array<T, 9> parameters = {
        k,         p,         v0,        normalization, fuseTime,
        threshold, ricochet0, ricochet1, nonAP};

    std::copy_n(reinterpret_cast<char *>(&parameters[0]),
                sizeof(T) * parameters.size(), &hashString[0]);
    return utility::base85Encode(hashString);
}

std::string generateHash(const double caliber, const double v0, const double cD,
                         const double mass, const double krupp,
                         const double normalization, const double fuseTime,
                         const double threshold, const double ricochet0,
                         const double ricochet1, const double nonAP) {
    double k = combinedAirDrag(cD, caliber, mass);
    double p = combinedPenetration(krupp, mass, caliber);
    return generateHash(k, p, v0, normalization, threshold, fuseTime, ricochet0,
                        ricochet1, nonAP);
}

std::string generateShellParamHash(
    const double caliber, const double v0, const double cD, const double mass,
    const double krupp, const double normalization, const double fuseTime,
    const double threshold, const double ricochet0, const double ricochet1,
    const double nonAP) {
    return generateHash(caliber, v0, cD, mass, krupp, normalization, fuseTime,
                        threshold, ricochet0, ricochet1, nonAP);
}

std::string generateHash(const shell &s) {
    return generateHash(s.caliber, s.v0, s.cD, s.mass, s.krupp, s.normalization,
                        s.fuseTime, s.threshold, s.ricochet0, s.ricochet1,
                        s.nonAP);
}

std::string generateShellHash(const shell &s) { return generateHash(s); }

}  // namespace wows_shell