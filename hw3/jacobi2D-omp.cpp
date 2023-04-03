// OpenMP hello world
// works with or without OPENMP
// $ g++ -fopenmp jacobi2D-omp.cpp && ./a.out
// $ g++ jacobi2D-omp.cpp $$ ./a.out
// depending on if you compile with or without openmp, the program will run

#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0; }
inline omp_int_t omp_get_num_threads() { return 1; }
#endif

#include <stdio.h>

#include <cmath>
#include <cstdlib>
#include <iostream>

void jacobi2D(double** u_new, double** u, const double** f, int N, double h,
              int num_iterations, int number_of_threads) {
    for (int k = 0; k < num_iterations; ++k) {
#pragma omp parallel for collapse(2) num_threads(number_of_threads)
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                double left = (i == 0) ? 0.0 : u[i - 1][j];
                double right = (i == N - 1) ? 0.0 : u[i + 1][j];
                double down = (j == 0) ? 0.0 : u[i][j - 1];
                double up = (j == N - 1) ? 0.0 : u[i][j + 1];

                u_new[i][j] =
                    0.25 * (h * h * f[i][j] + left + right + down + up);
            }
        }

        // Swap the pointers u_new and u for the next iteration
        std::swap(u_new, u);
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
    }  // number of threads

    double h = 1.0 / (N + 1);  // grid spacing

    double** u = new double*[N];
    double** u_new = new double*[N];
    double** f = new double*[N];

    for (int i = 0; i < N; ++i) {
        u[i] = new double[N];
        u_new[i] = new double[N];
        f[i] = new double[N];
    }

    // Initialize u and f
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            u[i][j] = 0;  // Initial values for u
            f[i][j] = 1;  // Set all elements of f to 1
        }
    }

#if defined(_OPENMP)
    double start_time = omp_get_wtime();
#endif

    // Call the Jacobi method with OpenMP
    jacobi2D(u_new, u, (const double**)f, N, h, num_iterations, num_threads);

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
        "Running time for %d iterations of jacobi 2d with size %d and %d "
        "threads is %f seconds\n",
        num_iterations, N, num_threads,
        end_time - start_time);  // Print the time (in seconds)
#endif

    // Deallocate memory
    for (int i = 0; i < N; ++i) {
        delete[] u[i];
        delete[] u_new[i];
        delete[] f[i];
    }
    delete[] u;
    delete[] u_new;
    delete[] f;
}
