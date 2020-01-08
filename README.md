# ShellCPP
Cpp header only library for calculating shell performance in World of Warships
## Compatibility
Should be compatible with most platforms - code is written to take advantage of multithreading and vectorization though 
## Current Functionality:
Calculate [For Various Distances]: 
### At Impact
angle of impact, impact velocity, raw belt/deck penetration, penetration adjusted for angle of impact
### Post-Penetration
shell detonation distances - some assumptiongs were made with regards to normalization changing shell direction
## Extensions
Python, WebAssembly
## Significant Issues:
Math library used has issues compiling with GCC / MSVC - some sort of inlining issue - should be fixed not test
## Disclaimer:
Apologies ahead of time if this repo appears to be made by an amateur (because I am one) - advice welcomed though
