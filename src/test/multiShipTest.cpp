#include <array>
#include <vector>

#include "../shellCPP.hpp"

int main() {
    std::vector<wows_shell::shell> ships;
    // ships.emplace_back(.460, 780, .292, 1460, 2574, 6, .033, 76, 45, 60, 0,
    //                   "Yamato");
    ships.emplace_back(.130, 870, .2857, 33.50, 1700, 10., .010, 22, 45, 60.0,
                       0, "Leningrad");
    ships.emplace_back(.203, 853, .3210, 118.0, 2846, 7.0, .033, 34, 60, 67.5,
                       0, "New Orleans");
    ships.emplace_back(.152, 950, .3210, 55.00, 2216, 8.5, .025, 25, 45, 60.0,
                       0, "Budyonny");
    ships.emplace_back(.220, 985, .2549, 176.0, 2590, 7.0, .033, 37, 45, 60.0,
                       0, "Moskva");
    ships.emplace_back(.150, 960, .3307, 45.50, 1862, 8.5, .025, 25, 45, 60.0,
                       0, "Nurnberg");
    ships.emplace_back(.283, 910, .3333, 300.0, 2282, 6.0, .010, 47, 45, 60.0,
                       0, "Graf Spee");
    ships.emplace_back(.152, 841, .3297, 50.8, 2609, 8.5, 0.005, 12, 60, 75.0,
                       0, "Edinburgh");
    ships.emplace_back(.152, 1000, .3256, 50, 2142, 8.5, 0.025, 25, 45, 60.0, 0,
                       "Duca d'Aosta");
    ships.emplace_back(.356, 792, .332, 680.4, 2604, 6, 0.033, 59, 45, 60.0, 0,
                       "Arizona");
    ships.emplace_back(.406, 701, .352, 1225, 2598, 6, .033, 68, 45, 60.0, 0,
                       "North Carolina");
    ships.emplace_back(.283, 890, .2827, 330, 2312, 6, 0.01, 47, 45, 60.0, 0,
                       "Scharnhorst");
    ships.emplace_back(.42, 800, .2994, 1220, 2415, 6, .033, 70, 45, 60.0, 0,
                       "Grosser Kurfurst 420");
    ships.emplace_back(.381, 731.5, .3379, 879, 2190, 6, 0.033, 64, 45, 60.0, 0,
                       "Hood");
    ships.emplace_back(.356, 757, .3142, 721, 2295, 6, .015, 59, 45, 60.0, 0,
                       "King George V");
    ships.emplace_back(.457, 762, .2897, 1506, 2485, 6, .033, 76, 45, 60.0, 0,
                       "Thunderer");
    ships.emplace_back(.33, 870, .2801, 560, 2428, 6, .033, 55, 45, 60.0, 0,
                       "Dunkerque");
    ships.emplace_back(.305, 762, .4595, 470.9, 2627, 6, 0.01, 51, 45, 60.0, 0,
                       "Oktyabrskaya Revolutsiya");
    ships.emplace_back(.32, 830, .4098, 525, 2600, 6, 0.033, 53, 45, 60.0, 0,
                       "Giulio Cesare");

    wows_shell::shellCalc calculator;
    // calculator.set_dt_min(.1);
    // calculator.set_precision(.01);
    std::array<double, 3> dists = {5000, 10000, 15000};
    for (wows_shell::shell &s : ships) {
        std::cout << "Ship Hash String: " << wows_shell::generateHash(s) << "\n";
        calculator
            .calculateImpact<false, wows_shell::numerical::forwardEuler>(s);
        std::cout << s.name << "\n";
        for (double dist : dists) {
            /*std::cout << s.interpolateDistanceImpact(
                             dist, wows_shell::impact::impactIndices::impactVelocity)
                      << " ";
            std::cout << s.interpolateDistanceImpact(
                             dist, wows_shell::impact::impactIndices::
                                       impactAngleHorizontalRadians)
                      << "; ";*/
            std::cout << s.interpolateDistanceImpact(
                             dist, wows_shell::impact::impactIndices::
                                       effectivePenetrationHorizontalNormalized)
                      << " ";
        }
        std::cout << "\n";
        
    }

    // ships.emplace_back(, , , , , , , , 45, 60.0,
    //                   0, "");
}