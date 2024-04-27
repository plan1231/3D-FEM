# FEM: 3D heat diffusion 

## Dependencies
This project assumes that [metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) and OpenMP
is already installed on your system. We use metis to partition input meshes into overlapping
subdomains, and OpenMP for paralellizing the program. On debian-based systems, run
```
apt install libomp-dev libmetis-dev
```
## Source code
All code lives under `/src`. 

`main.cpp` contains the bulk of the implementation,
from the Additive Schwarz Preconditioner to the simulation steps using implicit
Euler integration.

`domain_partitioning.cpp` contains the code for decomposing an input mesh (which
has to be a tetrahedron) into overlapping subdomains. We then build the restriction
and extension matrices for each subdomain.

## Compile

Compile this project using the standard cmake routine:
```
    mkdir build
    cd build
    cmake ..
    make
 ```
 the compiled binary should be called `3d_fem`. Use `3d_fem -h` to see the list of all options. 
 
 When running the program, your present working directory must be in `build`, as the path to the input files are hardcoded
 (see main.cpp:197). This might be fixed in the future.