# on Edision we will benchmark you against the default vendor-tuned BLAS. The compiler wrappers handle all the linking. If you wish to compare with other BLAS implementations, check the NERSC documentation.
# This makefile is intended for the GNU C compiler. To change compilers, you need to type something like: "module swap PrgEnv-pgi PrgEnv-gnu" See the NERSC documentation for available compilers.

#CC = gcc
CC = cc
OPT = -O3
OPT = -Ofast
#OPT = -march=corei7-avx
CFLAGS = $(OPT) -Wall -std=gnu99 -msse -msse2 -msse3 
LDFLAGS = -Wall 
# librt is needed for clock_gettime
LDLIBS = -lrt 
#LDLIBS = -framework Accelerate
targets = benchmark-naive benchmark-blocked benchmark-blas benchmark-final
objects = benchmark.o dgemm-naive.o dgemm-blocked.o dgemm-blas.o dgemm-final.o

.PHONY : default
default : all

.PHONY : all
all : clean $(targets)

benchmark-naive : benchmark.o dgemm-naive.o 
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-blocked : benchmark.o dgemm-blocked.o
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-blas : benchmark.o dgemm-blas.o
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-final: benchmark.o dgemm-final.o
	$(CC) -o $@ $^ $(LDLIBS)

%.o : %.c
	$(CC) -c $(CFLAGS) $<

.PHONY : clean
clean:
	rm -f $(targets) $(objects)
