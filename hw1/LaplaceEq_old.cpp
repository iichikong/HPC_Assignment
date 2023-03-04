
// problem3
// + g++ -O3 -std=c++11 LaplaceEq.cpp &&./a.out -n 10 -maxiter 1000

#include <math.h>
#include <stdio.h>

#include "utils.h"

int main(int argc, char** argv) {
    long N = read_option<long>("-n", argc, argv, "1");
    const long ite = read_option<long>("-maxiter", argc, argv, "5000");

    double h = 1.0 / (N + 1);
    double h2 = h * h;
    double* A = (double*)malloc(N * N * sizeof(double));  // initialization of A
    double* u = (double*)malloc(N * sizeof(double));      // initialication of u
    double* mem = (double*)malloc(N * sizeof(double));  // memory of previous u

    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            if (i == k) {
                A[i + k * N] = h2 * 2;
            } else if (i == k + 1 || i == k - 1) {
                A[i + k * N] = h2 * -1;
            } else {
                A[i + k * N] = 0;
            }
        }
    }

    double ini_res = sqrt(N);

    Timer t;

    // Jacobi method

    int k_jacobi;

    for (int i = 0; i < N; i++) {
        u[i] = 0;
        mem[i] = 0;
    }

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
                    mem[i] += A[i * N + j] * u[j];
                }
            }
            mem[i] = 1.0 / A[i * N + i] * (1.0 - mem[i]);
        }

        for (int i = 0; i < N; i++) {
            u[i] = mem[i];
            mem[i] = 0.0;
        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                mem[i] += A[i * N + j] * u[j];
            }
            res += (mem[i] - 1.0) * (mem[i] - 1.0);
        }
        res = sqrt(res);
        printf("%10d %10f\n", k, res);
        if (res * 10000 < ini_res) {
            k_jacobi = k;
            break;
        }
    }

    for (int i = 0; i < N; i++) {
        printf("%10f\n", u[i]);
    }

    double time_jacobi = t.toc();

    // Gauss-Seidel method

    for (int i = 0; i < N; i++) {
        u[i] = 0;
    }

    int k_gauss;

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
                    temp += A[i * N + j] * u[j];
                    // printf("%10d %10d %10f\n", i, j, A[i * N + j]);
                }
            }
            u[i] = 1.0 / A[i * N + i] * (1.0 - temp);
        }

        for (int i = 0; i < N; i++) {
            mem[i] = 0.0;
        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                mem[i] += A[i * N + j] * u[j];
            }
            res += (mem[i] - 1.0) * (mem[i] - 1.0);
        }
        res = sqrt(res);
        printf("%10d %10f\n", k, res);
        if (res * 10000 < ini_res) {
            k_gauss = k;
            break;
        }
    }

    double time_Gauss = t.toc();

    // runtime report
    printf("  run times for Jaccobi method: %10f\n", time_jacobi);
    printf("  number of iterations for Jaccobi method: %10ld\n", k_jacobi);

    printf("  run times for Gauss Seidel method: %10f\n", time_Gauss);
    printf("  number of iterations for Jaccobi method: %10ld\n", k_gauss);

    free(A);
    free(u);
    free(mem);

    return 0;
}