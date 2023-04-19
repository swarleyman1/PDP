CC = mpicc
CCFLAGS = -O3
LIBS = -lm

BINS = A2 test

all: $(BINS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)
