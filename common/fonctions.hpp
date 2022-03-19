#ifndef MATRICEVECTEUR_FONCTIONS_H
#define MATRICEVECTEUR_FONCTIONS_H

#include <chrono>

#define NOW std::chrono::system_clock::now()

void generate_vector(int n, int *vector, int zeros);
void matrix_vector(int n, int *matrix, int *v1, int *v2);

using Time = std::chrono::time_point <std::chrono::system_clock>;

#endif