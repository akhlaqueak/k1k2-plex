# Builds the serial and parallel programs: `make`, `make versions`, `make test`, `make clean`.
all:
	$(MAKE) -C serial
	$(MAKE) -C parallel

versions:
	$(MAKE) -C serial versions

test: all
	./tests/run_tests.sh

clean:
	$(MAKE) -C serial clean
	$(MAKE) -C parallel clean

.PHONY: all versions test clean
