cmake_minimum_required(VERSION 3.10)
project(Untangler3D CXX)

# --- Configuration ---

# Set C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# Find OpenCL
find_package(OpenCL REQUIRED)

# Find OpenMP
find_package(OpenMP)

set(GEOGRAM_ROOT "<YOUR GEOGRAM ROOT>")
set(GEOGRAM_BUILD "${GEOGRAM_ROOT}/build/Linux64-gcc-dynamic-Release")

# Common compilation flags (equivalent to your CXXFLAGS)
add_compile_options(-Wall)

# --- Targets ---
add_executable(untangler3d main.cpp untangler3d.cpp)
if(OpenMP_CXX_FOUND)
    target_link_libraries(untangler3d PUBLIC OpenMP::OpenMP_CXX)
endif()

target_include_directories(untangler3d PRIVATE
    "${GEOGRAM_ROOT}/src/lib"
    "${GEOGRAM_BUILD}/../include"
)

target_link_directories(untangler3d PRIVATE "${GEOGRAM_BUILD}/lib")

target_link_libraries(untangler3d PUBLIC OpenCL::OpenCL PRIVATE geogram )


target_compile_options(untangler3d PRIVATE 
    -O3 -march=native -ffast-math -funroll-loops -ftree-vectorize
)


add_custom_target(run
    COMMAND ./untangler3d ../domain_02.msh ../tangled.msh 
    DEPENDS untangler3d
    COMMENT "Running..."
    VERBATIM
)
