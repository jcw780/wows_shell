# ShellCPP
Cpp header only library for calculating shell performance in World of Warships
## Compatibility
Current version requires AVX2 supporting CPUs [most Intel/AMD architectures released post 2015 should work]
## Current Functionality:
Calculate [For Various Distances]: 
### At Impact
angle of impact, impact velocity, raw belt/deck penetration, penetration adjusted for angle of impact
### Post-Penetration
shell detonation distances
## Significant Issues:
Math library used has issues compiling with GCC / MSVC - some sort of inlining issue