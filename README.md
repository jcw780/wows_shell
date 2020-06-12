# ShellCPP
[![License](https://img.shields.io/github/license/jcw780/ShellCPP)](./LICENSE)
- C++ header only library for calculating shell performance in World of Warships
## Legal
- All copyrighted material provided by Wargaming.net are owned by Wargaming.net.
- All other material is available under the MIT License.
## Original Information / Source Code:
- Code: https://pastebin.com/1NEwkf7R
- Formulas: https://www.reddit.com/r/WorldOfWarships/comments/560yg2/wows_ballistic_model_penetration/
## Current Functionality:
### Shell Flight
Shell flight path
### At Impact:
Angle of Impact, Impact Velocity, Raw Belt/Deck Penetration, Penetration Adjusted for Angle of Impact and Normalization
Shell flight time (Real / In-game)
- Added ability to change trajectory computation method 
- Choices: - Forward Euler (original) - Runge-Kutta 2 and 4 - Adams-Bashforth 5
### Lateral Angles: 
- Definition: Angles where that represent the horizontal angling of a ship:
  + Bow-in: 90 degrees Full-broadside: 0 degrees
- Maximum Lateral Angle for Penetration - Minimum Lateral Angle for Fusing - Ricochet Lateral Angles
- Adjusts for angle of fall and vertical armor inclination
### Post-Penetration:
Shell detonation distance after penetration at various ranges, ship angling, and armor vertical inclinations
- Some assumptiongs were made with regards to normalization changing shell direction - testing is needed
- Added ability to modify the way the calculations are done 
### Fitting:
- Some limited capability to fit shells to real world data generating air drag and krupp values using gradient descent. 
## Compatibility:
- Requires C++17 supporting compiler
- Should be compatible with most platforms (though not extensively tested)
- Takes advantage of multithreading and vectorization (hardware supporting these should perform better if enabled)
## Extensions:
### Python 
- Written with Pybind11 - tested with Python 3.7.4 (Anaconda)
### WebAssembly 
- Uses Embind from Emscripten
- Could use threading - though has not been tested since browsers do not natively support wasm multithreading due to Spectre/Meltdown 
- Used in https://github.com/jcw780/wows_ballistics
## Future Goals:
- Makefiles
- Examples / Test Code Reorganization



