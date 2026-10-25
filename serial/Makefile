# Serial k1k2-plex. Override the compiler with e.g. `make CXX=clang++`.
CXX ?= g++
CXXFLAGS ?= -std=c++20 -O3
CC = $(CXX) -I../common/ $(CXXFLAGS)

all: k1k2plex

debug: CXXFLAGS += -g
debug: k1k2plex

# Ours: DTP + SOP (CTCP) + BRA + LAP
k1k2plex: kplex.cpp
	printf "#define ITERATIVE_PRUNE \n #define CTCP \n #define BRANCHING \n #define LOOKAHEAD" > version.h
	${CC}  kplex.cpp -o kplex

# Baselines (Base-0..3) and ablation variants of the paper
versions:
	printf "#define NAIVE" > version.h
	${CC}  kplex.cpp -o kplex-v0

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

	printf "#define TWO_HOP \n #define CTCP \n #define BRANCHING \n #define LOOKAHEAD" > version.h
	${CC}  kplex.cpp -o kplex-no-dtp

	printf "#define ITERATIVE_PRUNE \n #define CTCP \n #define LOOKAHEAD" > version.h
	${CC}  kplex.cpp -o kplex-no-branching

	printf "#define ITERATIVE_PRUNE \n #define CTCP \n #define BRANCHING" > version.h
	${CC}  kplex.cpp -o kplex-no-lookahead

	printf "#define ITERATIVE_PRUNE \n #define CTCP \n #define BRANCHING \n #define LOOKAHEAD" > version.h

clean:
	rm -f kplex kplex-v0 kplex-v1 kplex-v2 kplex-bk kplex-full kplex-no-ctcp kplex-no-dtp kplex-no-branching kplex-no-lookahead

.PHONY: all debug versions clean
