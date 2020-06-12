# Makefile for building the C language shared library for the l2p program

l2p.so: l2p.c
	R CMD SHLIB -o l2p.so -IR=1 -Wall -Wextra l2p.c rinterface.c l2pstats.c pwgenes.c utilfuncs.c a2a.c a2asupport.c

clean:
	-rm *.o *.so

