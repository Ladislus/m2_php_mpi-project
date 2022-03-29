#include "fonctions.hpp"

#define MATRIX_MAX 20

int main(int argc, char **argv) {

    // Init MPI
    int pid = -1, nprocs = -1;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    assert(nprocs == 1);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    assert(pid == 0);

    LOG(std::clog << "Master " << pid << "/" << nprocs - 1 << " started" << std::endl)

    // Read CLI arguments
    const unsigned int n = atoi(argv[1]);
    const unsigned int m = atoi(argv[2]);
    const unsigned int root = atoi(argv[4]);
    assert(root < nprocs);
    const unsigned int nslave = atoi(argv[5]);
    const std::string name = argv[6];
    const std::string slave_name = argv[7];

    LOG(
        if (pid == root) {
            std::clog << "\tn = " << n << std::endl;
            std::clog << "\tm = " << m << std::endl;
            std::clog << "\troot = " << root << std::endl;
            std::clog << "\tnslave = " << nslave << std::endl;
            std::clog << "\tname = '" << name << "'" << std::endl;
            std::clog << "\tslave_name = '" << slave_name << "'" << std::endl;
        }
    )

    // Define size constants
    const size_t nn = n * n;
    const size_t mn = m * n;

    // Spawn slaves & merge Comm
    LOG(std::clog << "[M" << pid << "] spawning " << nslave << " slaves" << std::endl)

    // Spawning slaves & merging Comm
    int intrapid = -1;
    MPI_Comm intercom = nullptr, intracom = nullptr;
    MPI_Comm_spawn(slave_name.c_str(), argv, nslave,
                   MPI_INFO_NULL, root, MPI_COMM_WORLD,
                   &intercom, MPI_ERRCODES_IGNORE);
    LOG(
        std::clog << "[M" << pid << "] spawning " << nslave << " slaves ✔" << std::endl;
        std::clog << "[M" << pid << "] merging comm" << std::endl
    )

    // FIXME: Infinite wait here, should match slave.cpp:39
    // Depending on which PC (and thus maybe which MPI implementation) I use, it might or might not work
    MPI_Intercomm_merge(intercom, 0, &intracom);
    MPI_Comm_rank(intracom, &intrapid);

    LOG(
        std::clog << "[M" << pid << "] merging comm ✔" << std::endl;
        std::clog << "[M" << pid << "] Generating matrix" << std::endl
    )

    // Initialize & broadcast matrix
    int *matrix = new int[nn];
    srand(time(nullptr));
    for (size_t i = 0; i < nn; i++) matrix[i] = rand() % MATRIX_MAX;

    LOG(
        std::clog << "[M" << pid << "] Generating matrix ✔" << std::endl;
//        for (size_t i = 0; i < n; ++i) {
//            std::clog << "\t" << i << ": [ ";
//            for (size_t j = 0; j < n; ++j) std::clog << matrix[i * n + j] << " ";
//            std::clog << "]" << std::endl;
//        }
        std::clog << "Broadcasting the matring to slaves" << std::endl
    )

    MPI_Bcast(matrix, nn, MPI_INT, root, intracom);

    LOG(
        std::clog << "Broadcasting the matring to slaves ✔" << std::endl;
        std::clog << "[M" << pid << "] genrating vectors" << std::endl
    )

    // initialize result and offset
    int offset = 0;
    int *results = new int[mn];

    // Initialize and generate vectors
    int *vectors = new int[mn];
    for (size_t i = 0; i < m; i++) generate_vector(n, vectors + (i * n), rand() % (n / 2));

    LOG(
        std::clog << "[M" << pid << "] genrating vectors ✔" << std::endl;
//        for (size_t i = 0; i < m; ++i) {
//            std::clog << "\tVector " << i <<" [ ";
//            for (size_t j = 0; j < n; ++j) std::clog << vectors[i * n + j] << " ";
//            std::clog << "]" << std::endl;
//        }
        std::clog << "[M" << pid << "] allocating windows" << std::endl
    )

    // Allocate windows
    MPI_Win vectors_win = nullptr, results_win = nullptr, offset_win = nullptr;
    MPI_Win_create(vectors, mn, sizeof(int), MPI_INFO_NULL, intracom, &vectors_win);
    MPI_Win_create(results, mn, sizeof(int), MPI_INFO_NULL, intracom, &results_win);
    MPI_Win_create(&offset, 1, sizeof(int), MPI_INFO_NULL, intracom, &offset_win);

    LOG(
        std::clog << "[M" << pid << "] allocating windows ✔" << std::endl;
        std::clog << "[M" << pid << "] waiting for first fence" << std::endl
    )

    // Fence to wait for windows initialization
    MPI_Win_fence(MPI_MODE_NOPRECEDE, vectors_win);

    LOG(
        std::clog << "[M" << pid << "] waiting for first fence ✔" << std::endl;
        std::clog << "[M" << pid << "] starting chrono" << std::endl
    )

    // Start chrono while slaves fetch & compute
    Time debut = NOW;

    LOG(
        std::clog << "[M" << pid << "] starting chrono ✔" << std::endl;
        std::clog << "[M" << pid << "] waiting for second fence" << std::endl
    )

    // Fence to wait for all vectors to be computed
    MPI_Win_fence(MPI_MODE_NOSUCCEED, results_win);

    LOG(std::clog << "[M" << pid << "] computing done" << std::endl)

    // Display computation time
    std::chrono::duration<double> elapsed_seconds = NOW - debut;
    std::cout << "Temps d'exécution: " << elapsed_seconds.count() << "s" << std::endl;

    // Open provided file
    std::fstream f(name, std::fstream::out);
    if (!f.is_open()) {
        std::cerr << "Error opening file '" << name << "'" << std::endl;
        std::cerr << strerror(errno) << std::endl;
        return EXIT_FAILURE;
    }

    // Write result into the file
    f << "Result:" << std::endl;
    for (size_t i = 0; i < m; i++) {
        f << "\tVector " << i << " [ ";
        for (size_t j = 0; j < n; j++) f << results[i * n + j] << " ";
        f << "]" << std::endl;
    }
    f.close();

//    LOG(
//        std::clog << "Result:" << std::endl;
//        for (size_t i = 0; i < m; i++) {
//            std::clog << "\tVector " << i << " [ ";
//            for (size_t j = 0; j < n; j++) std::clog << results[i * n + j] << " ";
//            std::clog << "]" << std::endl;
//        }
//    )

    // Free memory
    MPI_Comm_free(&intercom);
    MPI_Comm_free(&intracom);

    MPI_Win_free(&vectors_win);
    MPI_Win_free(&results_win);
    MPI_Win_free(&offset_win);

    MPI_Finalize();

    delete[] matrix;
    delete[] vectors;
    delete[] results;

    return EXIT_SUCCESS;
}