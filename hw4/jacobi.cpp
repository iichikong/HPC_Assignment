#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>

double compute_residual(double *lu, int lN, double invhsq){
    int i, j;
    double tmp, gres = 0.0, lres = 0.0;

    for (i = 1; i <= lN; i++){
        for (j = 1; j <= lN; j++) {
            tmp = ((4.0 * lu[i * (lN + 2) + j] - lu[(i - 1) * (lN + 2) + j] - lu[(i + 1) * (lN + 2) + j] - lu[i * (lN + 2) + j - 1] - lu[i * (lN + 2) + j + 1]) * invhsq - 1);
            lres += tmp * tmp;
        }
    }

    MPI_Allreduce(&lres, &gres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sqrt(gres);
}

int main(int argc, char * argv[]){
    int mpirank, i, j, p, N, lN, iter, max_iters;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    /* get name of host running MPI process */
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    printf("Rank %d/%d running on %s.\n", mpirank, p, processor_name);


    sscanf(argv[1], "%d", &N);
    sscanf(argv[2], "%d", &max_iters);

    int sqrt_p = (int)sqrt(p);

    if (N % sqrt_p != 0) {
        if (mpirank == 0) {
            printf("N: %d, sqrt_p: %d\n", N, sqrt_p);
            printf("Exiting. N must be a multiple of sqrt_p\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    lN = N / sqrt_p;

    MPI_Barrier(MPI_COMM_WORLD);
    double tt = MPI_Wtime();

    double * lu    = (double *) calloc(sizeof(double), (lN + 2) * (lN + 2));
    double * lunew = (double *) calloc(sizeof(double), (lN + 2) * (lN + 2));
    double * lutemp;

    double h = 1.0 / (N + 1);
    double hsq = h * h;
    double invhsq = 1./hsq;
    double gres, gres0, tol = 1e-5;

    gres0 = compute_residual(lu, lN, invhsq);
    gres = gres0;

    int row = mpirank / sqrt_p;
    int col = mpirank % sqrt_p;

    for (iter = 0; iter < max_iters && gres / gres0 > tol; iter++) {
        for (i = 1; i <= lN; i++) {
            for (j = 1; j <= lN; j++) {
                lunew[i * (lN + 2) + j]  = 0.25 * (hsq + lu[(i - 1) * (lN + 2) + j] + lu[(i + 1) * (lN + 2) + j] + lu[i * (lN + 2) + j - 1] + lu[i * (lN + 2) + j + 1]);
            }
        }

        // Communicate ghost values between neighbors
        int up = row > 0 ? mpirank - sqrt_p : MPI_PROC_NULL;
        int down = row < sqrt_p - 1 ? mpirank + sqrt_p : MPI_PROC_NULL;
        int left = col > 0 ? mpirank - 1 : MPI_PROC_NULL;
        int right = col < sqrt_p - 1 ? mpirank + 1 : MPI_PROC_NULL;

        // Send to up and receive from down
        MPI_Send(lunew + (lN + 2), lN, MPI_DOUBLE, up, 0, MPI_COMM_WORLD);
        MPI_Recv(lunew + (lN + 1) * (lN + 2), lN, MPI_DOUBLE, down, 0, MPI_COMM_WORLD, &status);

        // Send to down and receive from up
        MPI_Send(lunew + lN * (lN + 2), lN, MPI_DOUBLE, down, 1, MPI_COMM_WORLD);
        MPI_Recv(lunew, lN, MPI_DOUBLE, up, 1, MPI_COMM_WORLD, &status);

        // Send to left and receive from right
        for (i = 1; i <= lN; ++i) {
            MPI_Send(lunew + i * (lN + 2) + 1, 1, MPI_DOUBLE, left, 2, MPI_COMM_WORLD);
            MPI_Recv(lunew + i * (lN + 2) + lN + 1, 1, MPI_DOUBLE, right, 2, MPI_COMM_WORLD, &status);
        }

        // Send to right and receive from left
        for (i = 1; i <= lN; ++i) {
            MPI_Send(lunew + i * (lN + 2) + lN, 1, MPI_DOUBLE, right, 3, MPI_COMM_WORLD);
            MPI_Recv(lunew + i * (lN + 2), 1, MPI_DOUBLE, left, 3, MPI_COMM_WORLD, &status);
        }

        // Pointer swapping
        lutemp = lu; lu = lunew; lunew = lutemp;

        if (0 == (iter % 100)) {
            gres = compute_residual(lu, lN, invhsq);
            if (0 == mpirank) {
                printf("Iter %d: Residual: %g\n", iter, gres);
            }
        }
    }

    // Clean up
    free(lu);
    free(lunew);

    // Timing
    MPI_Barrier(MPI_COMM_WORLD);
    double elapsed = MPI_Wtime() - tt;
    if (0 == mpirank) {
        printf("Time elapsed is %f seconds.\n", elapsed);
    }
    MPI_Finalize();
    return 0;
}