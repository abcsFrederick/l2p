
CC=gcc
RM=rm -f
# fave ..
CFLAGS=-Wextra -Wall -O2 -funroll-loops -march=native -flto -pipe -minline-all-stringops 
LDFLAGS=-flto -lm

#debug
#CFLAGS=-g -Wall -Werror
#LDFLAGS=-lm

#CFLAGS=-Werror=implicit-function-declaration -Werror=implicit-fallthrough=3 -Werror=maybe-uninitialized -Werror=missing-field-initializers -Werror=incompatible-pointer-types -Werror=int-conversion -Werror=redundant-decls -Werror=parentheses -Wformat-nonliteral -Wformat-security -Wformat -Winit-self -Wmaybe-uninitialized -Wold-style-definition -Wredundant-decls -Wstrict-prototypes  -O2 -funroll-loops -march=native -flto -pipe -minline-all-stringops 

#debug
#CFLAGS=-g -fstack-protector -D_FORTIFY_SOURCE=2 -Wall
#LDFLAGS=-lm

#CFLAGS=-Wextra -Wall -Ofast -fomit-frame-pointer -march=native -flto -funroll-loops -pipe
#CFLAGS=-Wall -g -O2 -march=native -flto  -pipe
#LDFLAGS=-flto 

#CFLAGS=-Wall -O2
#LDFLAGS=

#CFLAGS=-g -pg -O3 -DL2PUSETHREADS=1 -Wall 
#LDFLAGS=-pg

#CFLAGS=-Ofast -fprofile-generate -fomit-frame-pointer -march=native -flto -funroll-loops
#LDFLAGS=-flto -fprofile-generate

#CFLAGS=-Ofast -fprofile-use -fomit-frame-pointer -march=native -flto -funroll-loops
#LDFLAGS=-flto -fprofile-use

LDLIBS=-lm -lpthread

SRCS=l2p.c pwgenes.c l2pstats.c utilfuncs.c  mitlicstats.c

OBJS=$(subst .c,.o,$(SRCS))

all: l2p

#
l2p: $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o l2p $(OBJS) $(LDLIBS)

l2p.o: small.h pathworks.h l2p.c utilfuncs.c 

pwgenes.o: small.h pathworks.h pwgenes.c

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) l2psmall


