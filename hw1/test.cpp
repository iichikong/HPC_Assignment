#include <stdio.h>

#include "utils.h"
double* MMult0(long m, long n, long k, double* a, double* b, double* c) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int p = 0; p < k; p++) {
                double A_ip = a[i + p * m];
                double B_pj = b[p + j * k];
                double C_ij = c[i + j * m];
                C_ij = C_ij + A_ip * B_pj;
                c[i + j * m] = C_ij;
            }
        }
    }
    return c;
}

int main(int argc, char** argv) {
    int p = 2;
    long m = p, n = p, k = p;
    double* a = (double*)malloc(m * k * sizeof(double));  // m x k
    double* b = (double*)malloc(k * n * sizeof(double));  // k x n
    double* c = (double*)malloc(m * n * sizeof(double));  // m x n

    // Initialize matrices
    for (long i = 0; i < m * k; i++) a[i] = 2;
    for (long i = 0; i < k * n; i++) b[i] = 2;
    for (long i = 0; i < m * n; i++) c[i] = 2;

    c = MMult0(m, n, k, a, b, c);
    for (int i = 0; i < n * n; i++) {
        printf("%10f", c[i]);
    }
    free(a);
    free(b);
    free(c);
    return 0;
}
