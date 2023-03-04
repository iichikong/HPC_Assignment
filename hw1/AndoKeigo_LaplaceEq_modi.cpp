// $ g++ -O3 -std=c++11 AndoKeigo_LaplaceEq_modi.cpp &&valgrind
// --tool=cachegrind cache-sim=yes ./a.out -n 10 -maxiter 5000

#include <math.h>
#include <stdio.h>

#include "utils.h"

// function for computing error
double error(long double* mem, long double* u, const long double a_ii,
             const long double a_ij, const long N) {
    long double res = 0.0;
    for (int i = 0; i < N; i++) {
        mem[i] = 0.0;
    }
    for (int i = 0; i < N; i++) {
        if (i == 0) {
            mem[i] = a_ii * u[i] + a_ij * u[i + 1];
        } else if (i == N - 1) {
            mem[i] += a_ij * u[i - 1] + a_ii * u[i];
        } else {
            mem[i] += a_ij * u[i - 1] + a_ii * u[i] + a_ij * u[i + 1];
        }
        res += (mem[i] - 1.0) * (mem[i] - 1.0);
    }
    res = sqrt(res);

    return res;
}

int main(int argc, char** argv) {
    long N = read_option<long>("-n", argc, argv, "1");
    const long ite = read_option<long>("-maxiter", argc, argv, "5000");

    long double h = 1.0 / (N + 1);
    long double* u =
        (long double*)malloc(N * sizeof(long double));  // initialication of u
    long double* mem = (long double*)malloc(
        N * sizeof(long double));  // temporary memory of previous u

    long double a_ii = h * h * 2;
    long double a_ij = h * h * -1;

    long double ini_res = sqrt(N);  // initial residual

    Timer t;

    ///////////////////
    // Jacobi method //
    ///////////////////

    int k_jacobi = ite;

    // initialize value of u and mem
    for (int i = 0; i < N; i++) {
        u[i] = 0;
        mem[i] = 0;
    }

    printf("\nJacobi method for N = %ld\n", N);
    printf("interation      u\n");
    printf("%10d %10Lf\n", 0, ini_res);

    t.tic();

    for (int k = 1; k < ite + 1; k++) {
        long double res = 0.0;

        for (int i = 0; i < N; i++) {
            mem[i] = 0.0;
            if (i == 0) {
                mem[i] = a_ij * u[i + 1];
            } else if (i == N - 1) {
                mem[i] += a_ij * u[i - 1];
            } else {
                mem[i] += a_ij * u[i - 1] + a_ij * u[i + 1];
            }

            mem[i] = 1.0 / a_ii * (1.0 - mem[i]);
        }
        for (int i = 0; i < N; i++) {
            u[i] = mem[i];
        }

        res = error(mem, u, a_ii, a_ij, N);

        printf("%10d %10Lf\n", k, res);

        if (res * 10000 < ini_res) {
            k_jacobi = k;
            break;
        }
    }

    double time_jacobi = t.toc();

    /////////////////////////
    // Gauss-Seidel method //
    /////////////////////////

    int k_gauss = ite;

    // initialize value of u
    for (int i = 0; i < N; i++) {
        u[i] = 0.0;
    }

    printf("\nGauss-Seidel method for N = %ld\n", N);
    printf("interation      u\n");

    printf("%10d %10Lf\n", 0, ini_res);

    t.tic();

    for (int k = 1; k < ite + 1; k++) {
        long double res = 0.0;

        for (int i = 0; i < N; i++) {
            long double temp = 0.0;
            if (i == 0) {
                temp = a_ij * u[i + 1];
            } else if (i == N - 1) {
                temp += a_ij * u[i - 1];
            } else {
                temp += a_ij * u[i - 1] + a_ij * u[i + 1];
            }
            u[i] = 1.0 / a_ii * (1.0 - temp);
        }

        res = error(mem, u, a_ii, a_ij, N);

        printf("%10d %10Lf\n", k, res);
        if (res * 10000 < ini_res) {
            k_gauss = k;
            break;
        }
    }

    double time_Gauss = t.toc();

    // runtime report
    printf("  run times for Jaccobi method: %10f\n", time_jacobi);
    printf("  the number of iterations for Jaccobi method: %10d\n", k_jacobi);

    printf("  run times for Gauss Seidel method: %10f\n", time_Gauss);
    printf("  the number of iterations for Gauss Seidel method: %10d\n",
           k_gauss);

    free(u);
    free(mem);

    return 0;
}