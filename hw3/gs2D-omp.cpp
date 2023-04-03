// OpenMP hello world
// works with or without OPENMP
// $ g++ -fopenmp gs2D-omp.cpp && ./a.out
// $ g++ gs2D-omp.cpp $$ ./a.out
// depending on if you compile with or without openmp, the program will run

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <stdio.h>

#include <cmath>
#include <cstdlib>
#include <iostream>

void gs2D(double** u, const double** f, int N, double h, int num_iterations,
          int num_threads) {
    for (int k = 0; k < num_iterations; ++k) {
        // Update red points

#pragma omp parallel for num_threads(num_threads)
        for (int i = 0; i < N; ++i) {
            for (int j = (i % 2); j < N; j += 2) {
                double left = (i > 0) ? u[i - 1][j] : 0;
                double right = (i < N - 1) ? u[i + 1][j] : 0;
                double up = (j > 0) ? u[i][j - 1] : 0;
                double down = (j < N - 1) ? u[i][j + 1] : 0;
                u[i][j] = (h * h * f[i][j] + left + right + up + down) / 4;
            }
        }

        // Update black points
#pragma omp parallel for num_threads(num_threads)
        for (int i = 0; i < N; ++i) {
            for (int j = ((i + 1) % 2); j < N; j += 2) {
                double left = (i > 0) ? u[i - 1][j] : 0;
                double right = (i < N - 1) ? u[i + 1][j] : 0;
                double up = (j > 0) ? u[i][j - 1] : 0;
                double down = (j < N - 1) ? u[i][j + 1] : 0;
                u[i][j] = (h * h * f[i][j] + left + right + up + down) / 4;
            }
        }
    }
}

int main(int argc, char* argv[]) {
    int N = 10;                 // Default size
    int num_iterations = 1000;  // Default number of iterations
    int num_threads = 8;

    // Parse command-line arguments
    if (argc > 1) {
        N = std::atoi(argv[1]);
    }
    if (argc > 2) {
        num_iterations = std::atoi(argv[2]);
    }  // size of the grid
    if (argc > 3) {
        num_threads = std::atoi(argv[3]);
    }  // size of the threads

    double h = 1.0 / (N + 1);  // grid spacing

    // Allocate memory for u and f
    double** u = new double*[N];
    double** f = new double*[N];
    for (int i = 0; i < N; ++i) {
        u[i] = new double[N];
        f[i] = new double[N];
    }

    // Initialize u and f
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            u[i][j] = 0;  // Initial values for u

            // Example function f(x, y) = sin(pi * x) * sin(pi * y)
            double x = (i + 1) * h;
            double y = (j + 1) * h;
            f[i][j] = 1;
        }
    }

#if defined(_OPENMP)
    double start_time = omp_get_wtime();
#endif

    // Call the Gauss-Seidel Red-Black function
    gs2D(u, (const double**)f, N, h, num_iterations, num_threads);

    if (N <= 10) {
        // Print the solution
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                std::cout << u[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

#if defined(_OPENMP)
    double end_time = omp_get_wtime();  // Record the end time

    printf(
        "Running time for %d iterations of gs 2d with size %d and %d threads "
        "is %f seconds\n",
        num_iterations, N, num_threads,
        end_time - start_time);  // Print the time (in seconds)
#endif

    // Deallocate memory for u and f
    for (int i = 0; i < N; ++i) {
        delete[] u[i];
        delete[] f[i];
    }
    delete[] u;
    delete[] f;

    return 0;
}