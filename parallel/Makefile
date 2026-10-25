# Parallel (OpenMP) k1k2-plex. Override e.g. `make CXX=clang++ OMPFLAGS="-Xpreprocessor -fopenmp -lomp"`.
CXX ?= g++
OMPFLAGS ?= -fopenmp
ifeq ($(shell uname -m),x86_64)
ARCHFLAGS ?= -msse4.2
endif
CXXFLAGS ?= -std=c++20 -O3 -Ofast
CC = $(CXX) -I../common/ $(CXXFLAGS) $(ARCHFLAGS) $(OMPFLAGS)

all: k1k2plex

debug: CXXFLAGS += -g
debug: k1k2plex

k1k2plex: kplex.cpp
	${CC}  kplex.cpp -o kplex

clean:
	rm -f kplex

.PHONY: all debug clean
