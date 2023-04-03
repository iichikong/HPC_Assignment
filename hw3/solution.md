## HPC
## Homework assignment 3
### Keigo Ando (ka2705)

1. 

- (a) For the first loop in `thread1` it takes `1+2+...+n/2=1/8n^2+n/4` and in `thread2`, it takes `n/2+(n/2+1)+...+n = 3/8n^2+1/4n` then the wating time is `1/4n^2`. For the second loop in `thread1` it takes `(n-1) + (n-2) + .... + n/2 = 3/8n^2 - 1/4n` and in `thread2` it takes `n/2 + (n/2 - 1) + ... + 0 = n^2/8 + n/4` in `thread2`. Then the waiting time for second loop is `1/4n^2-n/2`. Then total waiting time is `1/2n^2-1/2n`.

- (b) For the first loop in `thread1` it takes `1+3+...+(n-1) = n^2/4` and in `thread2` it takes `2+4+...+n = n^2/4+n/2` then the wating time is `n/2`. For the second loop in `thread1` it takes `(n-1) + (n-3) + ... + 1 = n^2/4` and in `thread2` it takes `(n-2) + (n-4) + ... + 0 = n^2/4 - n/2`. Then waiting time for the second loop is `n/2`. Then total waiting time is `n`.

- (c) It improves the speeding. For example, if `n=10`, for the second loop, `thread1` executes `f(9), f(6),f(5), f(2), f(1)` and `thread2` executes `f(8), f(7), f(4), f(3),f(0)`.  In this case, the waiting time is only `1`, which is much smaller than the case in (a) and (b).

- (d) It is not possible to completely eliminate waiting time by any means, since the total processing time may be an odd number. However, as we saw in the previous problem, it is possible to maximally reduce latency as much as possible by using dynamic or guided processing and by keeping the final chunk size small. In this case, the overhead may be the time required to process which process is assigned to which thread, so it is important to balance chunk size and other factors.

2.

The processor I used for this experiment is an 8-core AMD Ryzen 7 4700 with Radeon Graphics. Cpu MHz is 1996.193, and the cache size is 512 KB. And the program is executed on the WSL2 environment. 

The timing for different number of threads are as follwos:

|  number of threads  | timing  |
| ----:               |----:           |
|       1             |      0.356674     |
|       2             |      0.275255     |
|       3             |     0.214229     |
|       4             |     0.225793     |
|       6             |     0.288435     |
|       8             |       0.343615     |

The performance is not getting significally improved in my environment.


3.
My implementation can be seen in the source file `jacobi2D-omp.cpp` and `gs2D-omp.cpp`

I ran the experiment on the servers named `crunchy1.cims.nyu.edu`, which has Four AMD Opteron 6272 (2.1 GHz) (64 cores).

The timing of Jacobi and Gauss-Seidel method with 1000 iterations for different size and different number of threads can be seen as follows. As you can see, for the smaller size problem, paralllelizing does not make big difference(even it can be worse compared to single thread). That may be because of the overhead of setting up thread. However, for the larger problem, the timing is getting faster as the number of threads is getting larger 

|  N     |  number of threads  | Jacobi Method  | Gauss-Seidel Method | 
| ----:  | ----:               |----:           |----:                |
|  10    |       1             |      0.002555     |     0.002747              |
|  10    |       2             |      0.004079     |     0.005820              |
|  10    |       3             |      0.004622     |    0.008264           |
|  10    |       4             |      0.004668     |    0.008676            |
|  10    |       6             |      0.005869     |    0.010436            |
|  10    |       8             |      0.006893     |    0.012750            |
| 100   |       1              |  0.088943 |  0.037666  |
| 100   |       2              |  0.047124 |  0.025460  |
| 100   |       3              |  0.033919 |   0.021899  |
| 100   |       4              |  0.027149 |  0.018202  |
| 100   |       6              |  0.020749 |  0.019343  |
| 100   |       8              |  0.018911 |  0.022087  |
| 1000 | 1|9.264992| 6.958910|
| 1000 | 2|4.695436|3.354745|
| 1000 | 3|3.052067| 2.077986|
| 1000 | 4|2.274521|1.566572|
| 1000 | 6|1.518202|1.040622|
| 1000 | 8|1.152622|0.741681|
| 1000 | 16| 0.587426|0.397318|
| 1000 | 32| 0.393655|0.410206|
| 1000 | 64|  3.360101| 5.723483|

