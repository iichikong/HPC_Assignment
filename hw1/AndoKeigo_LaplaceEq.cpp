
// $ g++ -O3 -std=c++11 AndoKeigo_LaplaceEq.cpp &&valgrind --tool=cachegrind
// cache-sim=yes ./a.out -n 10 -maxiter 5000

#include <math.h>
#include <stdio.h>

#include "utils.h"

// function for computing error
double error(double* mem, double* u, const double a_ii, const double a_ij,
             const long N) {
    double res = 0.0;
    for (int i = 0; i < N; i++) {
        mem[i] = 0.0;
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                mem[i] += a_ii * u[j];
            }
            if (i == j + 1 || i == j - 1) {
                mem[i] += a_ij * u[j];
            }
        }
        res += (mem[i] - 1.0) * (mem[i] - 1.0);
    }
    res = sqrt(res);

    return res;
}

int main(int argc, char** argv) {
    long N = read_option<long>("-n", argc, argv, "1");
    const long ite = read_option<long>("-maxiter", argc, argv, "5000");

    double h = 1.0 / (N + 1);
    double* u = (double*)malloc(N * sizeof(double));  // initialication of u
    double* mem =
        (double*)malloc(N * sizeof(double));  // temporary memory of previous u

    double a_ii = h * h * 2;
    double a_ij = h * h * -1;

    double ini_res = sqrt(N);  // initial residual

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
    printf("%10d %10f\n", 0, ini_res);

    t.tic();

    for (int k = 1; k < ite + 1; k++) {
        double res = 0.0;

        for (int i = 0; i < N; i++) {
            mem[i] = 0.0;
            for (int j = 0; j < N; j++) {
                if (i == j) {
                    continue;
                } else {
                    if (i == j + 1 || i == j - 1) {
                        mem[i] += a_ij * u[j];
                    }
                }
            }
            mem[i] = 1.0 / a_ii * (1.0 - mem[i]);
        }
        for (int i = 0; i < N; i++) {
            u[i] = mem[i];
        }

        res = error(mem, u, a_ii, a_ij, N);

        printf("%10d %10f\n", k, res);

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
        u[i] = 0;
    }

    printf("\nGauss-Seidel method for N = %ld\n", N);
    printf("interation      u\n");

    printf("%10d %10f\n", 0, ini_res);

    t.tic();

    for (int k = 1; k < ite + 1; k++) {
        double res = 0.0;

        for (int i = 0; i < N; i++) {
            double temp = 0.0;
            for (int j = 0; j < N; j++) {
                if (i == j) {
                    continue;
                } else {
                    if (i == j + 1 || i == j - 1) {
                        temp += a_ij * u[j];
                    }
                }
            }
            u[i] = 1.0 / a_ii * (1.0 - temp);
        }

        res = error(mem, u, a_ii, a_ij, N);

        printf("%10d %10f\n", k, res);
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