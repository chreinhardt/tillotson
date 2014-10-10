gfortran -c icosahedron.f

cc -o ballic ballic.c tipsy.c icosahedron.o -lm


to use:

./ballic 40000 63.5 3e10 >tmp

