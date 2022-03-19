#include <iostream>
#include <fstream>
#include <chrono>
#include <mpi.h>
#include <cassert>

#include "fonctions.hpp"

#define VERBOSE 0

int main(int argc, char **argv) {

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    int pid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    assert(argc == 5);

    const int n = atoi(argv[1]);
    const int m = atoi(argv[2]);
    const int root = atoi(argv[3]);
    const std::string name = argv[4];

    const size_t nn = n * n;
    const size_t mn = m * n;
    const size_t data_size = sizeof(int);

    Time debut;

    int *matrix = new int[nn];
    int *vectors;

    std::fstream f;
    if (pid == root) {
        srand(time(nullptr));
        for (size_t i = 0; i < nn; i++) matrix[i] = rand() % 20;

#if VERBOSE
        // Print Matrix
        cout << "Matrix : " << endl;
        for (size_t i = 0; i < n; ++i) {
            cout << "\t";
            for (size_t j = 0; j < n; ++j) cout << matrix[i * n + j] << " ";
            cout << "" << endl;
        }
#endif

        vectors = new int[mn];
        for (size_t i = 0; i < m; i++) generate_vector(n, vectors + (i * n), rand() % (n / 2));

#if VERBOSE
        // Print vectors
        cout << "Vectors : " << endl;
        for (size_t i = 0; i < m; ++i) {
            cout << "\tVector " << i <<" [ ";
            for (size_t j = 0; j < n; ++j) cout << vectors[i * n + j] << " ";
            cout << "]" << endl;
        }
#endif
    }

    MPI_Bcast(matrix, nn, MPI_INT, root, MPI_COMM_WORLD);

    int sendcount[nprocs];
    int displs[nprocs];
    int counter = 0;
    const int qte = (m / nprocs);
    const int rst = m % nprocs;
    for (int i = 0; i < nprocs; ++i) {
        displs[i] = counter;
        int x = (rst > i) ? (qte + 1) * n : qte * n;
        sendcount[i] = x;
        counter += x;
    }
    const size_t sz = sendcount[pid] / n;
    int *batch = new int[sendcount[pid]];

    MPI_Scatterv(vectors, sendcount, displs, MPI_INT,
                 batch, sendcount[pid], MPI_INT,
                 root, MPI_COMM_WORLD);

#if VERBOSE
    // Print received vectors
    for (size_t i = 0; i < sz; ++i) {
        cout << "PID: " << pid << "; Vector " << i <<" [ ";
        for (size_t j = 0; j < n; ++j) cout << batch[i * n + j] << " ";
        cout << "]" << endl;
    }
#endif

    if (pid == root) debut = NOW;

    for (size_t i = 0; i < sz; ++i) {
        int result[n];
        matrix_vector(n, matrix, batch + (i * n), result);
        memcpy(batch + (i * n), result, n * data_size);
    }

    MPI_Gatherv(batch, sendcount[pid], MPI_INT,
                vectors, sendcount, displs, MPI_INT,
                root, MPI_COMM_WORLD);

    if (pid == root) {
        std::chrono::duration<double> elapsed_seconds = NOW - debut;
        std::cout << "Temps d'exécution : " << elapsed_seconds.count() << std::endl;
    }

    if (pid == root) {
        f.open(name, std::fstream::out);

        // Write result into the file
        f << "Result:" << std::endl;
        for (size_t i = 0; i < m; i++) {
            f << "\tVector " << i << " [ ";
            for (size_t j = 0; j < n; j++) f << vectors[i * n + j] << " ";
            f << "]" << std::endl;
        }
        f.close();

#if VERBOSE
        // Print result
        cout << "Result:" << endl;
        for (size_t i = 0; i < m; i++) {
            cout << "\tVector " << i << " [ ";
            for (size_t j = 0; j < n; j++) cout << vectors[i * n + j] << " ";
            cout << "]" << endl;
        }
#endif
    }

    MPI_Finalize();

    delete[] matrix;
    delete[] batch;
    if (pid == root) delete[] vectors;

    return 0;
}