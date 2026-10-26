# Tillotson EOS library (part of the EOSlib family).
#
# The library needs no external dependencies. Build it with -DTILL_NOEXIT to make failing
# evaluations return TILL_SENTINEL/TILL_FAIL instead of aborting (recommended for use in EOSlib;
# required by testtillotson). -DTILL_VERBOSE prints a message on every failed evaluation.

objects = tillotson.o tillinitlookup.o tillsplint.o

defs = -DTILL_NOEXIT
CFLAGS ?= -O3 -Wall $(defs)
LIBS ?= -lm

all: libtillotson.a testtillotson

libtillotson.a: $(objects)
	$(AR) rcs $@ $(objects)

$(objects): tillotson.h

# Invariant tests; the exit status is the number of failures.
testtillotson: testtillotson.o $(objects)
	$(CC) $(CFLAGS) -o $@ testtillotson.o $(objects) $(LIBS)

test: testtillotson
	./testtillotson

# Micro-benchmark of the hot-path functions.
bench: bench.o $(objects)
	$(CC) $(CFLAGS) -o $@ bench.o $(objects) $(LIBS)

# Any other single-file tool in this directory: make <name> builds <name>.c against the library.
%: %.o $(objects)
	$(CC) $(CFLAGS) -o $@ $< $(objects) $(LIBS)

clean:
	rm -f $(objects) testtillotson.o bench.o libtillotson.a testtillotson bench

.PHONY: all test clean