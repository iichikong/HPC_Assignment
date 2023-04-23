#include <mpi.h>
#include <iostream>
#include <cstdlib>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if (rank == 0) {
            std::cout << "Usage: mpirun -np <number_of_processes> ./int_ring_block <N_loops>\n";
        }
        MPI_Finalize();
        return 0;
    }

    int N_loops = atoi(argv[1]);
    int array_size = 2 * 1024 * 1024 / sizeof(int); // 2 MBytes of integers
    int* data = (int*) malloc(array_size * sizeof(int));
    if (rank == 0) {
        for (int i = 0; i < array_size; i++) {
            data[i] = i;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    for (int loop = 0; loop < N_loops; loop++) {
        if (rank == 0) {
            MPI_Send(data, array_size, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(data, array_size, MPI_INT, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(data, array_size, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(data, array_size, MPI_INT, (rank + 1) % size, 0, MPI_COMM_WORLD);
        }
    }

    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;

    if (rank == 0) {
        std::cout << "Total time: " << elapsed_time << " s" << std::endl;
        std::cout << "Estimated bandwidth: " << ((static_cast<double>(array_size) * sizeof(int) * N_loops * size) / elapsed_time) << " bytes/s" << std::endl;
    }

    free(data);
    MPI_Finalize();
    return 0;
}
