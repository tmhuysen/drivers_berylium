# Drivers for Berylium and Berylium 2-

## The CMake build step

In order to set up CMake for this project, externally build CMake:

```
mkdir build && cd build
cmake .. -DCMAKE_PREFIX_PATH=<GQCP_CMAKE_PREFIX_PATH> -DGQCP_INSTALL_PREFIX=<GQCP_INSTALL_PREFIX>
```

where `<GQCP_CMAKE_PREFIX_PATH>` is the `CMAKE_PREFIX_PATH` that you provided when installing the main library `gqcp` and `<GQCP_INSTALL_PREFIX>` is the path where `gqcp` is installed.


## Compiling the target executable

After successfully building the CMake system, use

```
make 
```

and if no errors occur, everything went correctly!
