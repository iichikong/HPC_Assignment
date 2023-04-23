#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if (rank == 0) {
            printf("Usage: mpirun -np <number_of_processes> ./int_ring <N_loops>\n" );
        }
        MPI_Finalize();
        return 0;
    }

    int N_loops = atoi(argv[1]);

    int value = 0;
    if (rank == 0) {
        value = rank;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    for (int loop = 0; loop < N_loops; loop++) {
        if (rank == 0) {
            MPI_Send(&value, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&value, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&value, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            value += rank;
            MPI_Send(&value, 1, MPI_INT, (rank + 1) % size, 0, MPI_COMM_WORLD);
        }
    }

    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;

    if (rank == 0) {
        printf("Result after %d loops: %d\n", N_loops, value);
        printf("Total time: %f s\n", elapsed_time);
        printf("Estimated latency: %f s\n", (elapsed_time / (N_loops * size)));
    }

    MPI_Finalize();
    return 0;
}
