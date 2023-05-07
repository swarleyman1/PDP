############################
# Makefile for lab 4
# Author: Malin Kallen 2019
# Modified: 2023-04-27
############################

CC = mpicc
CCFLAGS = -g -O3
LIBS = -lm

BINS = p_search

all: $(BINS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)
