#include <string>

#include "fonctions.hpp"

int main(int argc, char **argv) {

    /*
     * Depending on which PC I use, it might not work while using RMA, so use simple 'MPI_Init' instead
     * Error:
     *      The OSC pt2pt component does not support MPI_THREAD_MULTIPLE in this release.
     *      Workarounds are to run on a single node, or to use a system with an RDMA
     *      capable network such as Infiniband.
     */
    // Init MPI
    // int provided;
    // MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Init(&argc, &argv);

    int pid = -1, intraprid = -1, nprocs = -1;
    MPI_Comm intercom = nullptr, intracom = nullptr;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    assert(nprocs >= 1);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    assert(pid >= 0 && pid < nprocs);

    LOG(
        std::clog << "Slave " << pid << "/" << nprocs - 1 << " started" << std::endl;
        std::clog << "[S" << pid << "] getting intracom" << std::endl
    )
    // Get communicator for intra-process communication through merge
    MPI_Comm_get_parent(&intercom);
    MPI_Intercomm_merge(intercom, 1, &intracom);
    MPI_Comm_rank(intracom, &intraprid);
    assert(intraprid >= 0);
//    LOG(std::clog << "[S" << pid << "] getting intracom ✔" << std::endl)

    // Read CLI arguments
    const unsigned int n = atoi(argv[2]);
    const unsigned int m = atoi(argv[3]);
    const unsigned int batch_sz = atoi(argv[4]);
    const unsigned int root = atoi(argv[5]);
    assert(root < nprocs);

//    LOG(
//        if (pid == root) {
//            std::clog << "\tn = " << n << std::endl;
//            std::clog << "\tm = " << m << std::endl;
//            std::clog << "\tbatch_sz = " << batch_sz << std::endl;
//            std::clog << "\troot = " << root << std::endl;
//        }
//    )

    // Define size constant
    const size_t nn = n * n;

    // Allocate matrix memory & fetch from master
    LOG(std::clog << "[S" << pid << "] fetching matrix from master" << std::endl)
    int *matrix = new int[nn];
    MPI_Bcast(matrix, nn, MPI_INT, root, intracom);
    LOG(
//        std::clog << "[S" << pid << "] fetching matrix from master ✔" << std::endl;
        std::clog << "[S" << pid << "] allocating dull windows" << std::endl
    )

    // Allocate batch memory
    int *batch = new int[batch_sz * n];

    // Initialize dull windows (to match master initialization)
    MPI_Win vectors_win = nullptr, results_win = nullptr, offset_win = nullptr;
    MPI_Win_create(nullptr, 0, 1, MPI_INFO_NULL, intracom, &vectors_win);
    MPI_Win_create(nullptr, 0, 1, MPI_INFO_NULL, intracom, &results_win);
    MPI_Win_create(nullptr, 0, 1, MPI_INFO_NULL, intracom, &offset_win);

    LOG(
//        std::clog << "[S" << pid << "] allocating dull windows ✔" << std::endl;
        std::clog << "[S" << pid << "] waiting for first fence" << std::endl
    )

    // Fence to wait for windows initialization
    MPI_Win_fence(MPI_MODE_NOPRECEDE, vectors_win);

    LOG(
//        std::clog << "[S" << pid << "] waiting for first fence ✔" << std::endl;
        std::clog << "[S" << pid << "] starting work" << std::endl
    )

    int offset = -1, new_offset = -1;
    // Infinite loop (break on first condition when no more vectors to process)
    while (true) {
        LOG(std::clog << "[S" << pid << "] getting offset from master" << std::endl)
        // Get offset from master
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, root, 0, offset_win);
        MPI_Get(&offset, 1, MPI_INT, root, 0, 1, MPI_INT, offset_win);
        // If offset is -1, something went wrong with the previous MPI_Get, but MPI_SUCCESS was returned
        LOG(
//            std::clog << "[S" << pid << "] getting offset from master ✔" << std::endl;
            std::clog << "[S" << pid << "] offset = " << offset << ", batch_sz = " << batch_sz << std::endl
        )
        assert(offset >= 0);
        // Break if no more vectors to process
        if (new_offset >= m - 1 || offset >= m - 1) {
            LOG(std::clog << "[S" << pid << "] offset >= M, breaking" << std::endl)
            MPI_Win_unlock(root, offset_win);
            break;
        }
        // Get quantity of vectors to process (if not enough, get all remaining)
        const size_t sz = (offset + batch_sz > m) ? m - offset : batch_sz;
        // if sz > batch_sz, the received buffer will be overflown
        assert(sz <= batch_sz);
        // Compute the new vector offset for the other slaves
        new_offset = offset + sz;
        LOG(std::clog << "[S" << pid << "] putting new offset " << new_offset << " to master" << std::endl)
        // Update the offset on master
        MPI_Put(&new_offset, 1, MPI_INT, root, 0, 1, MPI_INT, offset_win);
        MPI_Win_unlock(root, offset_win);
        LOG(std::clog << "[S" << pid << "] putting new offset " << new_offset << " to master ✔" << std::endl)

        // Fetch the batch of vectors to process
        LOG(std::clog << "[S" << pid << "] getting batch from master" << std::endl)
        MPI_Win_lock(MPI_LOCK_SHARED, root, 0, vectors_win);
        MPI_Get(batch, sz * n, MPI_INT, root, offset * n, sz * n, MPI_INT, vectors_win);
        MPI_Win_unlock(root, vectors_win);

        LOG(
            std::clog << "[S" << pid << "] getting batch from master ✔" << std::endl;
//            std::clog << "Batch:" << std::endl;
//            for (std::size_t i = 0; i < sz; ++i) {
//                std::clog << "\t" << i << ": [ ";
//                for (std::size_t j = 0; j < n; ++j) std::clog << batch[i * n + j] << " ";
//                std::clog << "]" << std::endl;
//            }
            std::clog << "[S" << pid << "] processing current batch" << std::endl
        )

        // Process the batch
        for (size_t i = 0; i < sz; ++i) {
            // Temporary vector to store the result, prevent writing to the same vector we're reading
            int *result = new int[n];
            // Compute the dot product of the current vector with the matrix
            matrix_vector(n, matrix, batch + i * n, result);
            // Put the temporary result in the batch
            memcpy(batch + i * n, result, n * sizeof(int));
        }

        LOG(
//            std::clog << "[S" << pid << "] processing current batch ✔" << std::endl;
//            std::clog << "Batch result :" << std::endl;
//            for (std::size_t i = 0; i < sz; ++i) {
//                std::clog << "\tBatch " << i << ": [ ";
//                for (std::size_t j = 0; j < n; j++) std::clog << batch[i * n + j] << " ";
//                std::clog << "]" << std::endl;
//            }
            std::clog << "[S" << pid << "] putting results to master" << std::endl
        )

        // Put the result in the results window of the master
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, root, 0, results_win);
        MPI_Put(&batch, sz * n, MPI_INT, root, offset, sz * n, MPI_INT, results_win);
        MPI_Win_unlock(root, results_win);
        LOG(std::clog << "[S" << pid << "] putting results to master ✔" << std::endl)
    }
    LOG(std::clog << "[S" << pid << "] work ✔" << std::endl)

    // Fence to wait for all vectors to be computed
    MPI_Win_fence(MPI_MODE_NOSUCCEED, results_win);

    // Free memory
    MPI_Comm_free(&intercom);
    MPI_Comm_free(&intracom);

    MPI_Win_free(&vectors_win);
    MPI_Win_free(&results_win);
    MPI_Win_free(&offset_win);

    MPI_Finalize();
    delete[] matrix;
    delete[] batch;

    return EXIT_SUCCESS;
}