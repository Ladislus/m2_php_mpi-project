#include "fonctions.hpp"

int main(int argc, char **argv) {

	// Init MPI With OpenMP Support (Not used in this version)
	int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

	// Get the rank of the process
	int pid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	// Read arguments from CLI
	const int n = atoi(argv[1]);
    const int m = atoi(argv[2]);
    const int root = atoi(argv[3]);
    assert(root < nprocs);
    const std::string name = argv[4];

	// Define constants
	const size_t nn = n * n;
    const size_t mn = m * n;
    const size_t data_size = sizeof(int);

    Time debut;

	// Compute the number of vectors per process
    int count[nprocs];
    int displs[nprocs];
    const int qte = m / nprocs;
    const int rst = m % nprocs;
    int cpt = 0;
    for (int i = 0; i < nprocs; ++i) {
        displs[i] = cpt;
        const int x = (rst > i) ? (qte + 1) * n : qte * n;
        count[i] = x;
        cpt += x;
    }
    const size_t sz = count[pid] / n;

    int *matrix = new int[nn];
    int *batch = new int[count[pid]];
    int *vectors;

	// Create the dynamic window
    long data_address = -1;
    MPI_Win window = nullptr;
    assert(MPI_Win_create_dynamic(MPI_INFO_NULL, MPI_COMM_WORLD, &window) == MPI_SUCCESS);

	// Root process attach the data to the window
    if (pid == root) {
        assert(MPI_Win_attach(window, matrix, nn * data_size) == MPI_SUCCESS);
        MPI_Get_address(matrix, &data_address);
        assert(data_address >= 0);
    }
	// Broadcast the data address to all processes for further use
    MPI_Bcast(&data_address, 1, MPI_LONG, root, MPI_COMM_WORLD);

	// Root process generates the matrix & vectors
	std::fstream f;
    if (pid == root) {
		// Init the matrix
        srand(time(nullptr));
        for (size_t i = 0; i < nn; i++) matrix[i] = rand() % 20;

        LOG(
            // Print Matrix
            std::clog << "Matrix : " << std::endl;
            for (size_t i = 0; i < n; ++i) {
                std::clog << "\t";
                for (size_t j = 0; j < n; ++j) std::clog << matrix[i * n + j] << " ";
                std::clog << std::endl;
            }
        )

		// Init the vectors
        vectors = new int[mn];
        for (size_t i = 0; i < m; i++) generate_vector(n, vectors + (i * n), rand() % (n / 2));

        LOG(
            // Print vectors
            std::clog << "Vectors : " << std::endl;
            for (size_t i = 0; i < m; ++i) {
                std::clog << "\tVector " << i <<" [ ";
                for (size_t j = 0; j < n; ++j) std::clog << vectors[i * n + j] << " ";
                std::clog << "]" << std::endl;
            }
        )
    }

	// Fence to wait for master to finish initialization
    MPI_Win_fence(MPI_MODE_NOPRECEDE, window);

    if (pid != root) {
        // Fetch the matrix
        MPI_Win_lock(MPI_LOCK_SHARED, root, 0, window);
        MPI_Get(matrix, nn, MPI_INT, root, data_address, nn, MPI_INT, window);
        MPI_Win_unlock(root, window);
    }

	// Fence to wait for threads to finish reading the matrix
    MPI_Win_fence(0, window);

	// Root process detach from matrix & attach to vectors
    if (pid == root) {
        MPI_Win_detach(window, matrix);
        assert(MPI_Win_attach(window, vectors, mn * data_size) == MPI_SUCCESS);
        MPI_Get_address(vectors, &data_address);
        assert(data_address >= 0);
    }
	// Broadcast the data address to all processes for further use
    MPI_Bcast(&data_address, 1, MPI_LONG, root, MPI_COMM_WORLD);

    if (pid != root) {
        // Fetch the vectors
        MPI_Win_lock(MPI_LOCK_SHARED, root, 0, window);
        MPI_Get(batch, count[pid], MPI_INT, root, data_address + displs[pid] * data_size, count[pid], MPI_INT, window);
        MPI_Win_unlock(root, window);
    } else
        // Copy vectors into batch
        memcpy(batch, vectors + displs[pid] * data_size, count[pid] * data_size);

    LOG(
        // Print received vectors
        for (size_t i = 0; i < sz; ++i) {
            std::clog << "PID: " << pid << "; Vector " << i <<" [ ";
            for (size_t j = 0; j < n; ++j) std::clog << batch[i * n + j] << " ";
            std::clog << "]" << std::endl;
        }
    )

	// Fence to wait for threads to finish reading the vectors
    MPI_Win_fence(0, window);

	// Save the starting time
	if (pid == root) debut = NOW;

	// Compute the dot product of the vectors and the matrix
	for (size_t i = 0; i < sz; ++i) {
        int result[n];
        matrix_vector(n, matrix, batch + (i * n), result);
        memcpy(batch + (i * n), result, n * data_size);
    }

	// Fence to wait for threads to finish computing the dot product
    MPI_Win_fence(0, window);

	// Threads send the result to the master
    if (pid != root) {
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, root, 0, window);
        MPI_Put(batch, count[pid], MPI_INT, root, data_address + displs[pid] * data_size, count[pid], MPI_INT, window);
        MPI_Win_unlock(root, window);
    } else
        memcpy(vectors + displs[pid] * data_size, batch, count[pid] * data_size);

	// Fence to wait for threads to finish writing the result to master
    MPI_Win_fence(0, window);

	// Root process output process time
	if (pid == root) {
        std::chrono::duration<double> elapsed_seconds = NOW - debut;
        std::cout << "Time: " << elapsed_seconds.count() << "s" << std::endl;
    }

	// Fence to wait for master to finish outputing the time
	// Most likely useless
    MPI_Win_fence(MPI_MODE_NOSUCCEED, window);

	// Root write the result to file
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

        LOG(
            // Print result
            std::clog << "Result:" << std::endl;
            for (size_t i = 0; i < m; i++) {
                std::clog << "\tVector " << i << " [ ";
                for (size_t j = 0; j < n; j++) std::clog << vectors[i * n + j] << " ";
                std::clog << "]" << std::endl;
            }
        )
    }

	// Finalize MPI & free memory
    MPI_Win_free(&window);
    MPI_Finalize();

    delete[] matrix;
    delete[] batch;
    if (pid == root) delete[] vectors;

    return 0;
}