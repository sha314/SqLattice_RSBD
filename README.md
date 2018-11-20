# SqLattice_RSBD
Random Sequential Ballistic Deposition on Square Lattice for studying long range phase transition

This repository was created only to reference my work for my thesis on phase transition

### How to compile
This program is CMake based and requires CMake version >= 3.
Also gnu make, g++, gcc is required.
The general building method of a CMake program is applied here.
The program must be build in release mode (default mode in general) for better performance.

### How to run the program
In the executable directory run the executable with two additional command line argument
first one is the length of the lattice and second one is the ensemble size.
The program will generate data frile for l={0,1,2} for each run but it has a clever method 
for not using the same file name for different run 
(it uses time to generate distinct data file name)
