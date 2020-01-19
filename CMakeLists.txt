# top-level CMakeLists.txt

cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(gqcp-link LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)

# Find gqcp and link the executable
set(CMAKE_MODULE_PATH "${GQCP_PREFIX_PATH}/cmake/;${GQCP_INSTALL_PREFIX}/cmake")
find_package(gqcp REQUIRED)

add_executable(berylium berylium.cpp)
target_link_libraries(berylium PUBLIC gqcp::gqcp)


# cmake .. -DCMAKE_PREFIX_PATH=/miniconda3 -DGQCP_INSTALL_PREFIX=/miniconda3 