
objects = table.o tillotson.o interpol/coeff.o interpol/interpol.o interpol/brent.o

CFLAGS ?= -O3

table: $(objects)
	cc -o table $(objects) -lm

all: table

clean:
	rm $(objects)

