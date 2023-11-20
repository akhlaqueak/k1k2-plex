CC=g++ -I../common/ -std=c++20 -fopenmp -O3

all: k1k2plex

debug: CC+= -g
debug: k1k2plex

k1k2plex: kplex.cpp
	git pull
	${CC}  kplex.cpp -o kplex

