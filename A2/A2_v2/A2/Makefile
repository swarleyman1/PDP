############################
# Makefile for lab 2
# Author: Malin Kallen 2019
############################

CC = mpicc
CCFLAGS = -g -O3
LIBS = -lm

BINS = matmul

all: $(BINS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)
