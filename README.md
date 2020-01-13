# ShellCPP
- C++ header only library for calculating shell performance in World of Warships
## Original Source:
- Code: https://pastebin.com/1NEwkf7R
- Formulas: https://www.reddit.com/r/WorldOfWarships/comments/560yg2/wows_ballistic_model_penetration/
## Current Functionality:
### Shell Flight
Shell flight path
### At Impact:
Angle of Impact, Impact Velocity, Raw Belt/Deck Penetration, Penetration Adjusted for Angle of Impact and Normalization
Shell flight time (Real / In-game)
### Post-Penetration:
Shell detonation distance after penetration at various ranges, ship angling, and armor vertical inclinations
- Some assumptiongs were made with regards to normalization changing shell direction - testing is needed
- Added ability to modify the way the calculations are done 
## Compatibility:
### Software:
C++17 
### Hardware: 
- Should be compatible with most platforms 
- Code does take advantage of multithreading and vectorization (platforms supporting these should perform better)
## Extensions:
### Python 
- Written with Pybind11 - tested with Python 3.7.4 (Anaconda)
### WebAssembly 
- Written with Embind - threading does not working since atomics do not work yet
- Used in https://jcw780.github.io/overpen_calculator/ 
## Future Goals:
- Makefiles
- Examples / Test Code Reorganization




