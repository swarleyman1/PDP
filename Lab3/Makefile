CC = mpicc
CCFLAGS = -g -O3
LIBS = -lm

BINS = array_sect array_test cart_shift_NSEW sendrecv_replace

all: $(BINS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)
