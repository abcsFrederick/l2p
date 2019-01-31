# Makefile for building the C language shared library for the l2p program

l2p.so: l2p.c
	R CMD SHLIB -o l2p.so -IR=1  l2p.c

clean:
	-rm *.o *.so

