
objects = tillotson.o interpol/coeff.o interpol/interpol.o interpol/brent.o

CFLAGS ?= -O3

default:
	@echo "Please specify which tool you want to make."

all:
	default

table: table.o $(objects)
	cc -o table table.o $(objects) -lm

cold: cold.o $(objects)
	cc -o cold cold.o $(objects) -lm

coldlookup: coldlookup.o $(objects)
	cc -o coldlookup coldlookup.o $(objects) -lm
clean:
	rm $(objects)

