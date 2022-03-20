#ifndef MATRICEVECTEUR_FONCTIONS_HPP
#define MATRICEVECTEUR_FONCTIONS_HPP

#include <iostream>
#include <chrono>
#include <fstream>
#include <cassert>
#include <mpi.h>

#include "config.h"


// Macro to include debug logs at compile time
#define LOG(code) if (VERBOSE) { code; }

/*
 * Macro to automatically add (at compile time) a test on MPI functions to keep code (somewhat) readable,
 * As they all return an int corresponding to the possible error code
 * MPI_SUCCESS is 0
 * 55 was chosen arbitrarily
 * (Code is wrapped inside bracket to prevent variable conflict/redefinition)
 * /!\ This is not a good practice, but it's a simple one
 * /!\ DO NOT ADD A SEMICOLON AFTER THE CODE INSIDE THE MACRO !
 */
/* Became useless as MPI check by default
#define MPI_ABORT_CODE 55
#define CHECK(code) {                                                                                                                                                       \
                        int __ret = (code);                                                                                                                                 \
                        if (__ret != MPI_SUCCESS) {                                                                                                                         \
                            std::cerr << "MPI code '" << #code << "' (" << __FILE__ << "::" << __LINE__ << ") errored with code " << __ret << ", aborting !" << std::endl;  \
                            MPI_Abort(MPI_COMM_WORLD, MPI_ABORT_CODE); }                                                                                                    \
                    }
*/

#define NOW std::chrono::system_clock::now()

void generate_vector(size_t size, int *vector, size_t zeros);
void matrix_vector(size_t size, const int *matrix, const int *vector, int *output);
void matrix_vector_omp(size_t size, const int *matrix, const int *vector, int *output);

using Time = std::chrono::time_point <std::chrono::system_clock>;

#endif