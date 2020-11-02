# ShellCPP
[![License](https://img.shields.io/github/license/jcw780/ShellCPP)](./LICENSE)
![GitHub commit activity](https://img.shields.io/github/commit-activity/y/jcw780/ShellCPP?style=plastic)
[![Discord](https://discordapp.com/api/guilds/731224331136532531/widget.png)](https://discord.gg/fpDB9y5)
- C++ header only library for calculating shell performance in World of Warships
- Experimentally verified to be nearly identical to the World of Warships ballistics model
## Legal
- All copyrighted material provided by Wargaming.net are owned by Wargaming.net.
- All other material is available under the MIT License.
## Original Information & Source Code:
Original Inspiration:
- Code: https://pastebin.com/1NEwkf7R
- Formulas: https://www.reddit.com/r/WorldOfWarships/comments/560yg2/wows_ballistic_model_penetration <br/>

Revisions:
- Code: https://pastebin.com/GXUt7BMJ
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
- Should be compatible with most platforms 
- Ability to utilize multithreading and vectorization for maximum performance
## Extensions:
### Python 
- Written with Pybind11
- Tested with Python 3.7.4 (Anaconda), 3.8.5
### WebAssembly 
- Uses Embind from Emscripten
- Used in https://github.com/jcw780/wows_ballistics
## Future Goals:
- Wiki / Tutorial
- Threading on WebAssembly module when brower support becomes commonplace
- Updating penetration formula once more refined values are acquired



