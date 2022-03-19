CMAKE_MINIMUM_REQUIRED(VERSION 3.13)
PROJECT("MPI Project"
        VERSION 1.0.0
        LANGUAGES CXX C)

# Set C++ standard (version) and add compile flags for warnings
SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_C_STANDARD 11)
ADD_COMPILE_OPTIONS(-O3 -Wall -Wextra -pedantic)

# Find common libraries sources and add the folder to the include path²
FILE(GLOB common "common/*.cpp")
INCLUDE_DIRECTORIES(common)

# Find MPI & OpenMPI
FIND_PACKAGE(OpenMP REQUIRED)
FIND_PACKAGE(MPI REQUIRED)
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")

# Vanilla version (MPI 1: Scatterv/Gatherv)
SET(VANILLA_EXECUTABLE Vanilla)
FILE(GLOB vanilla_sources "vanilla/*.cpp")
ADD_EXECUTABLE(${VANILLA_EXECUTABLE} ${vanilla_sources} ${common})
TARGET_LINK_LIBRARIES(${VANILLA_EXECUTABLE} MPI::MPI_CXX)

# Vanilla version with OpenMP
SET(OPENMP_EXECUTABLE OMP)
FILE(GLOB omp_sources "omp/*.cpp")
ADD_EXECUTABLE(${OPENMP_EXECUTABLE} ${omp_sources} ${common})
TARGET_LINK_LIBRARIES(${OPENMP_EXECUTABLE} MPI::MPI_CXX OpenMP::OpenMP_CXX)