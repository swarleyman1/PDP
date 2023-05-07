### Makefile for A3 ###


CC = mpicc
CCFLAGS = -g -O3
LIBS = -lm

BINS = quicksort

all: $(BINS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)
