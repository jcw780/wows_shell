# wows_shell
[![License](https://img.shields.io/github/license/jcw780/ShellCPP)](./LICENSE)
![GitHub commit activity](https://img.shields.io/github/commit-activity/y/jcw780/ShellCPP?style=plastic)
[![Discord](https://discordapp.com/api/guilds/731224331136532531/widget.png)](https://discord.gg/fpDB9y5)
- C++ header only library for calculating shell performance in World of Warships
- Experimentally verified to be nearly identical to the World of Warships ballistics model
## Legal
- All copyrighted material provided by Wargaming.net are owned by Wargaming.net.
- All other material is available under the MIT License.
## Original Model Information & Source Code:
Original Model:
- [Code](https://pastebin.com/1NEwkf7R)
- [Post-Penetration Formulas](https://www.reddit.com/r/WorldOfWarships/comments/560yg2/wows_ballistic_model_penetration) <br/>

Revised Model:
- [Code](https://pastebin.com/GXUt7BMJ)
## Features:
### Shell Flight
Computes shell flight path.
### At Impact:
Outputs:
  - Angle of Impact
  - Impact Velocity 
  - Raw, Belt/Deck, Normalization Adjusted Penetration 
  - Shell flight time (Real / In-game) <br/>

Ability to change trajectory computation method
  - Forward Euler (default - also in-game method) 
  - Runge-Kutta 2 and 4 
  - Adams-Bashforth 5
### Lateral Angles: 
Computes lateral angles where penetration, AP fuzing, and ricochets occur while adjusting for impact angle and vertical armor inclination. <br/>

Definition: Angles where that represent the horizontal angling of a ship:
  - Bow-in: 90 degrees 
  - Full-broadside: 0 degrees <br/>

Outputs:
  - Maximum Lateral Angle for Penetration 
  - Minimum Lateral Angle for Fusing 
  - Ricochet Lateral Angles <br/>

### Post-Penetration:
Shell detonation distance after penetration while adjusting for ship angling, and vertical armor inclinations. <br/>

Ability to modify the way the calculations are done 
- Enable or Disable Normalization changing direction
- Linear estimation or full air drag modeling
### Fitting:
Capability to fit shells to real world data using gradient descent. 
- Air Drag Coefficient
- Krupp
## Compatibility:
- Requires C++17 supporting compiler
- Supports multithreading and vectorization for improved performance
## Extensions:
### Python 
- Requires Pybind11
- Tested with Python 3.7.4 (Anaconda), 3.8.5
### WebAssembly 
- Requires Emscripten
- Used in https://github.com/jcw780/wows_ballistics
## Future Goals:
- Wiki / Tutorial
- Packages, more refined build tools
- Will continue to update the model when new information is acquired



