CC=g++ -I../common/ -std=c++20 -O3

all: k1k2plex

debug: CC+= -g
debug: k1k2plex

k1k2plex: kplex.cpp
	printf "#define ITERATIVE_PRUNE \n #define CTCP \n #define BRANCHING \n #define LOOKAHEAD" > version.h
	${CC}  kplex.cpp -o kplex

versions:
	printf "#define NAIVE" > version.h
	${CC}  kplex.cpp -o kplex-v0

	printf "#define NAIVE \n #define CTCP" > version.h
	${CC}  kplex.cpp -o kplex-v1

	printf "#define NAIVE \n #define CTCP" > version.h
	${CC}  kplex.cpp -o kplex-v1
	
	printf "#define TWO_HOP \n #define CTCP" > version.h
	${CC}  kplex.cpp -o kplex-v2
	
	printf "#define ITERATIVE_PRUNE \n #define CTCP" > version.h
	${CC}  kplex.cpp -o kplex-bk

	printf "#define ITERATIVE_PRUNE \n #define CTCP \n #define BRANCHING \n #define LOOKAHEAD" > version.h
	${CC}  kplex.cpp -o kplex-full

	printf "#define ITERATIVE_PRUNE \n #define BRANCHING \n #define LOOKAHEAD" > version.h
	${CC}  kplex.cpp -o kplex-no-ctcp

	printf "#define CTCP \n #define BRANCHING \n #define LOOKAHEAD" > version.h
	${CC}  kplex.cpp -o kplex-no-iterative

	printf "#define ITERATIVE_PRUNE \n #define CTCP \n #define LOOKAHEAD" > version.h
	${CC}  kplex.cpp -o kplex-no-branching

	printf "#define ITERATIVE_PRUNE \n #define CTCP \n #define BRANCHING" > version.h
	${CC}  kplex.cpp -o kplex-no-lookahead