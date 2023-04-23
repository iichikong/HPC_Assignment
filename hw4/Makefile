CC = mpicc
CFLAGS = -std=c11 -O3 -march=native -fopenmp
LIBS = -lm
TARGETS = int_ring int_ring_block jacobi pingpong

all : $(TARGETS)

int_ring: int_ring.c
	$(CC) $(CFLAGS) $< -o $@

int_ring_block: int_ring_block.c
	$(CC) $(CFLAGS) $< -o $@

jacobi: jacobi.c
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

pingpong: pingpong.c
	$(CC) $(CFLAGS) $< -o $@

clean:
		-$(RM) $(TARGETS) *~

.PHONY: all, clean