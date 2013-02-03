
default: ParaSparse

ParaSparse: ParaSparse.c ParaSparse.h PS_runner.c PS_runner.h
	mpicc -g -Wall -o ParaSparse ParaSparse.c PS_runner.c

clean:
	rm -f ParaSparse

