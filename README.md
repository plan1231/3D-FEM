# FEM: 3D heat diffusion 

## Dependencies
This project assumes that [metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) and OpenMP
are already installed on your system. We use metis to partition input meshes into overlapping
subdomains, and OpenMP for paralellizing the program. On debian-based systems, run
```
apt install libomp-dev libmetis-dev
```
## Source code
All code lives under `src`. 

`main.cpp` contains the bulk of the implementation,
from the Additive Schwarz Preconditioner to the simulation steps using implicit
Euler integration.

`domain_partitioning.cpp` contains the code for decomposing an input mesh (which
has to be a tetrahedron) into overlapping subdomains. We then build the restriction
and extension matrices for each subdomain.

## Compiling
For release mode with all optimizations
```
cmake --preset=gcc -DCMAKE_BUILD_TYPE=Release -B release
cd release
make
 ```

For debugging:
```
cmake  -DCMAKE_BUILD_TYPE=Debug -B build
cd build
make
```

## Tetgen errors
The version of Tetgen that ships with libigl seems to contain a bug. When compiling for the first time in Debug mode,
your compiler might throw this error:
```
In file included from /usr/include/c++/12/cassert:44,
                 from /root/3D-FEM/build/_deps/libigl-src/include/igl/copyleft/tetgen/mesh_to_tetgenio.cpp:14,
                 from /root/3D-FEM/build/_deps/libigl-src/include/igl/copyleft/tetgen/mesh_to_tetgenio.h:53,
                 from /root/3D-FEM/build/_deps/libigl-src/include/igl/copyleft/tetgen/tetrahedralize.cpp:9,
                 from /root/3D-FEM/build/_deps/libigl-src/include/igl/copyleft/tetgen/tetrahedralize.h:124,
                 from /root/3D-FEM/src/main.cpp:7:
/root/3D-FEM/build/_deps/libigl-src/include/igl/copyleft/tetgen/tetgenio_to_tetmesh.cpp: In function ‘bool igl::copyleft::tetgen::tetgenio_to_tetmesh(const tetgenio&, Eigen::PlainObjectBase<Derived>&, Eigen::PlainObjectBase<DerivedF>&, Eigen::PlainObjectBase<DerivedT>&, Eigen::PlainObjectBase<DerivedF>&, Eigen::PlainObjectBase<DerivedFTC>&, Eigen::PlainObjectBase<DerivedFN>&, Eigen::PlainObjectBase<DerivedTV>&, Eigen::PlainObjectBase<DerivedTT>&, int&)’:
/root/3D-FEM/build/_deps/libigl-src/include/igl/copyleft/tetgen/tetgenio_to_tetmesh.cpp:86:38: error: ‘const class tetgenio’ has no member named ‘numberofpointmarkers’; did you mean ‘numberofpointmtrs’?
   86 |     assert(out.numberofpoints == out.numberofpointmarkers);
      |                                      ^~~~~~~~~~~~~~~~~~~~
make[2]: *** [CMakeFiles/3d_fem.dir/build.make:90: CMakeFiles/3d_fem.dir/src/main.cpp.o] Error 1
make[1]: *** [CMakeFiles/Makefile2:428: CMakeFiles/3d_fem.dir/all] Error 2
make: *** [Makefile:156: all] Error 2
```
To fix this, go to the file that contains this error and delete this line.


 the compiled binary should be called `3d_fem`. Use `3d_fem -h` to see the list of all options. 
 
 When running the program, your present working directory must be in `build`, as the path to the input files are hardcoded
 (see main.cpp:197). This might be fixed in the future.