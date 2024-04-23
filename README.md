# FEM: 3D heat diffusion 

## Dependencies
This project assumes that [metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
is already installed on your system. We need it to partition input meshes into overlapping
subdomains. On a mac, simply run
```
brew install metis
```


## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `example` binary.

## Run

From within the `build` directory just issue:

    ./example

A glfw app should launch displaying a 3D cube.
