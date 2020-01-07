# ShellCPP
Cpp header only library for calculating shell performance in World of Warships
## Compatibility
~~Current version requires AVX2 supporting CPUs [most Intel/AMD architectures released post 2015 should work]
- Working on a version that does not use intrinsics - may alleviate compile issues with GCC and MSVC though uncertain about performance~~
- No intrinsic option added
## Current Functionality:
Calculate [For Various Distances]: 
### At Impact
angle of impact, impact velocity, raw belt/deck penetration, penetration adjusted for angle of impact
### Post-Penetration
shell detonation distances - some assumptiongs were made with regards to normalization changing shell direction
## Significant Issues:
Math library used has issues compiling with GCC / MSVC - some sort of inlining issue
## Disclaimer:
Apologies ahead of time if this repo appears to be made by an amateur (because I am one) - advice welcomed though
