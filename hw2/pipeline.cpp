// g++ -g pipeline.cpp -o pipeline && valgrind
// --leak-check=full ./pipeline
// g++ pipeline.cpp -o pipeline && ./pipeline

#include <stdio.h>
#include <stdlib.h>

#include "utils.h"

int main(int argc, char **argv) {
    Timer t;
    // long n = read_option<long>("-n", argc, argv);
    const long n = 400000000;

    double *a = (double *)malloc(n * sizeof(double));
    double *b = (double *)malloc(n * sizeof(double));
    double s = 0.0;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    double sum4 = 0.0;
    double temp1 = 0.0;
    double temp2 = 0.0;
    double temp3 = 0.0;
    double temp4 = 0.0;

    // Initialize matrices
    for (long i = 0; i < n; i++) a[i] = drand48();
    for (long i = 0; i < n; i++) b[i] = drand48();

    t.tic();
    for (long i = 0; i < n; i++) {
        s += a[i] * b[i];
    }
    printf("expe1: %f ms  %f\n", t.toc() * 1000, s);

    s = 0.0;
    t.tic();
    for (long i = 0; i < n / 2 - 1; i++) {
        sum1 += a[2 * i] * b[2 * i];
        sum2 += a[2 * i + 1] * b[2 * i + 1];
    }
    s = sum1 + sum2;
    printf("expe2: %f ms  %f\n", t.toc() * 1000, s);

    s = 0.0;
    sum1 = 0.0;
    sum2 = 0.0;
    t.tic();
    for (long i = 0; i < n / 2 - 1; i++) {
        sum1 += *(a + 0) * *(b + 0);
        sum2 += *(a + 1) * *(b + 1);
        a += 2;
        b += 2;
    }
    s = sum1 + sum2;
    printf("expe3: %f ms  %f\n", t.toc() * 1000, s);

    a -= (n - 2);
    b -= (n - 2);

    s = 0.0;
    sum1 = 0.0;
    sum2 = 0.0;
    t.tic();
    for (long i = 0; i < n / 2 - 1; i++) {
        temp1 = *(a + 0) * *(b + 0);
        temp2 = *(a + 1) * *(b + 1);

        sum1 += temp1;
        sum2 += temp2;

        a += 2;
        b += 2;
    }
    s = sum1 + sum2;
    printf("expe4: %f ms  %f\n", t.toc() * 1000, s);

    a -= (n - 2);
    b -= (n - 2);

    s = 0.0;
    sum1 = 0.0;
    sum2 = 0.0;
    t.tic();
    for (long i = 0; i < n / 2 - 1; i++) {
        sum1 += temp1;
        temp1 = *(a + 0) * *(b + 0);

        sum2 += temp2;
        temp2 = *(a + 1) * *(b + 1);

        a += 2;
        b += 2;
    }
    s = sum1 + sum2;
    printf("expe5: %f ms  %f\n", t.toc() * 1000, s);

    a -= (n - 2);
    b -= (n - 2);

    s = 0.0;
    sum1 = 0.0;
    sum2 = 0.0;
    t.tic();
    for (long i = 0; i < n / 4 - 1; i++) {
        sum1 += a[4 * i] * b[4 * i];
        sum2 += a[4 * i + 1] * b[4 * i + 1];
        sum3 += a[4 * i + 2] * b[4 * i + 2];
        sum4 += a[4 * i + 3] * b[4 * i + 3];
    }
    s = sum1 + sum2 + sum3 + sum4;
    printf("expe2: %f ms  %f\n", t.toc() * 1000, s);

    s = 0.0;
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    t.tic();
    for (long i = 0; i < n / 4 - 1; i++) {
        sum1 += *(a + 0) * *(b + 0);
        sum2 += *(a + 1) * *(b + 1);
        sum3 += *(a + 2) * *(b + 2);
        sum4 += *(a + 3) * *(b + 3);

        a += 4;
        b += 4;
    }
    s = sum1 + sum2 + sum3 + sum4;
    printf("expe3: %f ms  %f\n", t.toc() * 1000, s);

    a -= (n - 4);
    b -= (n - 4);

    s = 0.0;
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    t.tic();
    for (long i = 0; i < n / 4 - 1; i++) {
        temp1 = *(a + 0) * *(b + 0);
        temp2 = *(a + 1) * *(b + 1);
        temp3 = *(a + 2) * *(b + 2);
        temp4 = *(a + 3) * *(b + 3);

        sum1 += temp1;
        sum2 += temp2;
        sum3 += temp3;
        sum4 += temp4;

        a += 4;
        b += 4;
    }
    s = sum1 + sum2 + sum3 + sum4;
    printf("expe4: %f ms  %f\n", t.toc() * 1000, s);

    a -= (n - 4);
    b -= (n - 4);

    s = 0.0;
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    t.tic();
    for (long i = 0; i < n / 4 - 1; i++) {
        sum1 += temp1;
        temp1 = *(a + 0) * *(b + 0);

        sum2 += temp2;
        temp2 = *(a + 1) * *(b + 1);

        sum3 += temp3;
        temp3 = *(a + 2) * *(b + 2);

        sum4 += temp4;
        temp4 = *(a + 3) * *(b + 3);

        a += 4;
        b += 4;
    }
    s = sum1 + sum2 + sum3 + sum4;
    printf("expe5: %f ms  %f\n", t.toc() * 1000, s);

    a -= (n - 4);
    b -= (n - 4);

    free(a);
    free(b);
}