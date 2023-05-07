### Makefile for project ###


CC = mpicc
CCFLAGS = -g -O3
LIBS = -lm



mc: monte_carlo_SSA.c
	$(CC) $(CCFLAGS) monte_carlo_SSA.c -o mc $(LIBS)

clean:
	$(RM) mc