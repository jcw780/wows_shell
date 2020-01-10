# ShellCPP
- C++ header only library for calculating shell performance in World of Warships
## Original code for ballistic calculations and post penetration velocity formula comes from: 
- Code: https://pastebin.com/1NEwkf7R
- Post Penetration: https://www.reddit.com/r/WorldOfWarships/comments/560yg2/wows_ballistic_model_penetration/
## Current Functionality:
### Shell Flight
Shell flight path
### At Impact:
Angle of Impact, Impact Velocity, Raw Belt/Deck Penetration, Penetration Adjusted for Angle of Impact and Normalization
Shell flight time (Real / In-game)
### Post-Penetration:
Shell detonation distance after penetration at various ranges, ship angling, and armor vertical inclinations
- some assumptiongs were made with regards to normalization changing shell direction - testing is needed
## Compatibility:
Should be compatible with most platforms
Code is written to take advantage of multithreading and vectorization so platforms with such functionality (Multicore CPUs SSE/AVX support - basically every relatively recent desktop CPU) should perform better
## Extensions:
### Python 
- written with Pybind11 - tested with Python 3.7.4 (Anaconda)
### WebAssembly 
- written with Embind - threading does not working since atomics do not work yet
- used in https://jcw780.github.io/overpen_calculator/ 
## Future Goals:
- Makefiles
- Examples / Test Code Reorganization




