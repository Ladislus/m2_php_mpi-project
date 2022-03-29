# MPI Project
Project made during the second year of Master's Degree in Computer Science (IMIS) at the University of Orléans for the PHP module (High performance computing).

## Members
 - WALCAK Ladislas
 - QUETIER Thomas

## Requirements
 - MPI: We used OpenMPI version 4.0.2
 - OMP (For the OMP version): We used libomp-dev 0.3.6-1

## Launch
To launch the project, you can use the custom targets created in the CMakeLists.txt file (`run_XXX`, ie: `run_vanilla`, `run_omp`, ...). 
Arguments and number of processes are controlled by the following variables: 
 - **PROC_QUANTITY**: Define the number of processes to launch (-np option of MPI).
 - **DEFAULT_ARGUMENTS**: Define the default arguments for all targets except Master/Slave.
 - **MS_ARGUMENTS**: Define the arguments for Master/Slave targets (-np is ignored and set explicitly set 1, as multiple masters are useless).

## Codes
The project contains 5 versions:
 - [Vanilla](vanilla.cpp): The vanilla version of the project, using MPI 1 (Scatterv/Gatherv).
 - [OMP](omp.cpp): The vanilla version of the project, using MPI 1 (Scatterv/Gatherv) extended with OMP.
 - [RMA](rma.cpp): The RMA version of the project, using MPI windows.
 - [RMA_OMP](rma_omp.cpp): The RMA version of the project, using MPI windows extended with OMP.
 - [Master](master_slave/master.cpp) / [Slave](master_slave/slave.cpp): The MPI 3 version of the project, using MPI Master/Slave technique (MPI_Comm_spawn) with RMA windows.

The code contained in the [common](common) folder is shared between all versions.  
We didn't make a Master/Slave version using MPI 1 because we didn't thing it would be useful.  
We couldn't make a Master/Slave version extended with OMP, as it isn't possible to use the current version of OpenMP in a Master/Slave version (MPI_THREAD_MULTIPLE, required to use OpenMP, provoque an error with the current version).

## Ptimirev
This section is dedicated to the use of the project on the University's computing bay.

Ptimirev CMake isn't at least 3.11, so the project is not usable with CMake (We used CMake 3.11 import targets, so we cannot downgrade it). 
To use it, you need to use the [bench.sh](bench/bench.sh) script to compile and launch the benchmarks (you must be inside the [bench/](bench) folder for it to work). 
Because we used a [config.h.in](common/config.h.in) file, which is not generated without CMake, you must copy/paste the content of the config file inside the required files (fonctions.hpp), 
and remove the imports. You can either do it manually, or use the code on the`cmake` git branch, which is already set up to be used on ptimirev.