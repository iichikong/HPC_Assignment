#include <math.h>
#include <omp.h>
#include <stdio.h>

#include <algorithm>

// Scan A array and write result into prefix_sum array;
// use long data type to avoid overflow
void scan_seq(long* prefix_sum, const long* A, long n) {
    if (n == 0) return;
    prefix_sum[0] = 0;
    for (long i = 1; i < n; i++) {
        prefix_sum[i] = prefix_sum[i - 1] + A[i - 1];
    }
}

void scan_omp(long* prefix_sum, const long* A, long n) {
    int p = 8;  // omp_get_max_threads();
    int chunk_size = (n + p - 1) / p;

    // Step 1: Split the input array into p parts and perform local scans

#pragma omp parallel num_threads(p)
    {
        int t = omp_get_thread_num();

        int start = t * chunk_size;
        int end = std::min(start + chunk_size, (int)n);
        // printf("thread %d of %d: start %d\n", t + 1, p, start);

        long local_sum = 0;
        for (int i = start; i < end; ++i) {
            local_sum += A[i];
            prefix_sum[i + 1] = local_sum;
        }
        // #pragma pragma omp parallel for num_threads(p)
        // for (long i = 1; i < n; i++) {
        //     prefix_sum[i] = prefix_sum[i - 1] + A[i - 1];
        // }
    }

    // Step 2: Compute the correction values
    long* partial_sums = (long*)malloc(p * sizeof(long));
    partial_sums[0] = 0;
    for (int i = 1; i < p; ++i) {
        partial_sums[i] =
            partial_sums[i - 1] + prefix_sum[i * ((n + p - 1) / p)];
    }

// Step 3: Update the partial sums with the corrections
#pragma omp parallel num_threads(p)
    {
        int t = omp_get_thread_num();
        if (t > 0) {
            int start = t * chunk_size;
            int end = std::min(start + chunk_size, (int)n);
            long correction = partial_sums[t];
            for (int i = start; i < end; ++i) {
                prefix_sum[i + 1] += correction;
            }
        }
    }
    free(partial_sums);
}

int main() {
    long N = 100000000;

    long* A = (long*)malloc(N * sizeof(long));
    long* B0 = (long*)malloc(N * sizeof(long));
    long* B1 = (long*)malloc(N * sizeof(long));
    for (long i = 0; i < N; i++) A[i] = rand();
    for (long i = 0; i < N; i++) B1[i] = 0;

    double tt = omp_get_wtime();
    scan_seq(B0, A, N);
    printf("sequential-scan = %fs\n", omp_get_wtime() - tt);

    tt = omp_get_wtime();
    scan_omp(B1, A, N);
    printf("parallel-scan   = %fs\n", omp_get_wtime() - tt);

    long err = 0;
    for (long i = 0; i < N; i++) err = std::max(err, std::abs(B0[i] - B1[i]));
    printf("error = %ld\n", err);

    free(A);
    free(B0);
    free(B1);
    return 0;
}
