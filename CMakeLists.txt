# top-level CMakeLists.txt

cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(gqcp-link LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)

# Find gqcp and link the executable
set(CMAKE_MODULE_PATH "${GQCP_PREFIX_PATH}/cmake/;${GQCP_INSTALL_PREFIX}/cmake")
find_package(gqcp REQUIRED)

add_executable(berylium2 berylium2.cpp)
add_executable(berylium4 berylium4.cpp)
add_executable(berylium berylium0.cpp)
add_executable(prio prio.cpp)
target_link_libraries(berylium PUBLIC GQCP::gqcp)
target_link_libraries(berylium2 PUBLIC GQCP::gqcp)
target_link_libraries(berylium4 PUBLIC GQCP::gqcp)

# cmake .. -DCMAKE_PREFIX_PATH=~/miniconda3 -DGQCP_INSTALL_PREFIX=~/miniconda3 
