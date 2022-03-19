#include <iostream>
#include "fonctions.hpp"

void generate_vector(int n, int *vector, int zeros) {
    for (size_t i = 0; i < zeros; i++) vector[i] = 0;
    for (size_t i = zeros; i < n; i++) vector[i] = rand() % 20;
}

void matrix_vector(int n, int *matrix, int *v1, int *v2) {
    int ptr = 0;
    while (v1[ptr] == 0) ++ptr;
    for (size_t i = 0; i < n; i++) {
        v2[i] = 0;
        for (size_t j = ptr; j < n; j++) v2[i] += v1[j] * matrix[i * n + j];
    }
}