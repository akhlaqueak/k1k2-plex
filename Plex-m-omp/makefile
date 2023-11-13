FLAG = -O3 -std=c++14 -msse4.2 -Ofast -fopenmp

#FLAG = -g --std=c++14 -O2 
ifdef PROMPT
	FLAG := $(FLAG) -DPROMPT
endif

ifdef OPT
	FLAG := $(FLAG) -Ofast
endif
#-fcilkplus

all: 
	g++ $(FLAG) listPlex.cpp -o listPlex
	g++ toBin.cpp -o toBin
clean:
	rm listPlex
