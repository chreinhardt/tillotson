objects = tillotson.o tillinitlookup.o tillsplint.o interpol/brent.o nr/nrcubicspline.o nr/nrutil.o

exe = table pressureoldnew lookup lookup_cold testu1 testspline testsplint testnewsplint testsplint2 testsplinerho testsplintrho testsplinev testsplintv testcubicintrho testlookupucold testudrho testudv testgrid testpolyv printderiv printpress pressneg testisintable testisbelowcoldcurve testrhomin testoutofbounds testsolvebc calcisentrope testrhoptemp calcpressure testdirectintegration testpoverrho2 testeospressure testtillpressure testtillrhopu tillpressrhotemp

#defs = -DTILL_PRESS_NP -DTILL_OUTPUT_ALL_WARNINGS -DTILL_PRESS_MELOSH
defs = -DTILL_PRESS_NP -DTILL_OUTPUT_ALL_WARNINGS

CFLAGS ?= -O3 $(defs)

default:
	@echo "Please specify which tool you want to make."

all:
	default

# Debug the function tillLookupU() by comparing the results of the lookup to a
# direct calculation using tillCalcU() which simply solves the ODE using a RK4
# integrator.
table: table.o $(objects)
	cc -o table table.o $(objects) -lm

# Compare the results of tillPressureSoundold() that contains the old code we used in
# Gasoline and tillPressureSound()
pressureoldnew: pressureoldnew.o $(objects)
	cc -o pressureoldnew pressureoldnew.o $(objects) -lm

#
# Generate a lookup table for a given material.
#
lookup: lookup.o $(objects)
	cc -o lookup lookup.o $(objects) -lm

#
# Make a lookup table for the cold curve.
#
lookup_cold: lookup_cold.o $(objects)
	cc -o lookup_cold lookup_cold.o $(objects) -lm

testu1: testu1.o $(objects)
	cc -o testu1 testu1.o $(objects) -lm

testspline: testspline.o $(objects)
	cc -o testspline testspline.o $(objects) -lm

#
# Code for debugging the interpolation function tillCubicInt().
#
testsplint: testsplint.o $(objects)
	cc -o testsplint testsplint.o $(objects) -lm

#
# Code for debugging the new logarithmic lookup table.
#
testsplintlogrho: testsplintlogrho.o $(objects)
	cc -o testsplintlogrho testsplintlogrho.o $(objects) -lm

#
# Code for debugging tillLookupU using the new logarithmic lookup table.
#
testlookupulogrho: testlookupulogrho.o $(objects)
	cc -o testlookupulogrho testlookupulogrho.o $(objects) -lm

#
# Code for debugging tillCalcU using log(rho) as integration variable.
#
testcalcu: testcalcu.o $(objects)
	cc -o testcalcu testcalcu.o $(objects) -lm

#
# Pretty much the same but it compares the old with the new interpolator.
#
testnewsplint: testnewsplint.o $(objects)
	cc -o testnewsplint testnewsplint.o $(objects) -lm

testsplint2: testsplint2.o $(objects)
	cc -o testsplint2 testsplint2.o $(objects) -lm

testsplinerho: testsplinerho.o $(objects)
	cc -o testsplinerho testsplinerho.o $(objects) -lm

testsplintrho: testsplintrho.o $(objects)
	cc -o testsplintrho testsplintrho.o $(objects) -lm

testsplinev: testsplinev.o $(objects)
	cc -o testsplinev testsplinev.o $(objects) -lm

testsplintv: testsplintv.o $(objects)
	cc -o testsplintv testsplintv.o $(objects) -lm

testcubicintrho: testcubicintrho.o $(objects)
	cc -o testcubicintrho testcubicintrho.o $(objects) -lm

testlookupucold: testlookupucold.o $(objects)
	cc -o testlookupucold testlookupucold.o $(objects) -lm

testudrho: testudrho.o $(objects)
	cc -o testudrho testudrho.o $(objects) -lm

testudv: testudv.o $(objects)
	cc -o testudv testudv.o $(objects) -lm

testgrid: testgrid.o $(objects)
	cc -o testgrid testgrid.o $(objects) -lm

testpolyv: testpolyv.o $(objects)
	cc -o testpolyv testpolyv.o $(objects) -lm

printderiv: printderiv.o $(objects)
	cc -o printderiv printderiv.o $(objects) -lm

printpress: printpress.o $(objects)
	cc -o printpress printpress.o $(objects) -lm

pressneg: pressneg.o $(objects)
	cc -o pressneg pressneg.o $(objects) -lm

testisintable: testisintable.o $(objects)
	cc -o testisintable testisintable.o $(objects) -lm

testisbelowcoldcurve: testisbelowcoldcurve.o $(objects)
	cc -o testisbelowcoldcurve testisbelowcoldcurve.o $(objects) -lm

testrhomin: testrhomin.o $(objects)
	cc -o testrhomin testrhomin.o $(objects) -lm

testoutofbounds: testoutofbounds.o $(objects)
	cc -o testoutofbounds testoutofbounds.o $(objects) -lm

testsolvebc: testsolvebc.o $(objects)
	cc -o testsolvebc testsolvebc.o $(objects) -lm

calcisentrope: calcisentrope.o $(objects)
	cc -o calcisentrope calcisentrope.o $(objects) -lm

testrhoptemp: testrhoptemp.o $(objects)
	cc -o testrhoptemp testrhoptemp.o $(objects) -lm

calcpressure: calcpressure.o $(objects)
	cc -o calcpressure calcpressure.o $(objects) -lm

# Code to test the direct integration tillCalcU().
testdirectintegration: testdirectintegration.o $(objects)
	cc -o testdirectintegration testdirectintegration.o $(objects) -lm

#
# Print PoverRho2 for small rho to see, if the expression diverges.
#
testpoverrho2: testpoverrho2.o $(objects)
	cc -o testpoverrho2 testpoverrho2.o $(objects) -lm

#
# Test the function eosPressureSound.
#
testeospressure: testeospressure.o $(objects)
	cc -o testeospressure testeospressure.o $(objects) -lm

#
# Test the function eosdPdrho.
#
testeosdpdrho: testeosdpdrho.o $(objects)
	cc -o testeosdpdrho testeosdpdrho.o $(objects) -lm

#
# Test the function eosdPdu.
#
testeosdpdu: testeosdpdu.o $(objects)
	cc -o testeosdpdu testeosdpdu.o $(objects) -lm

#
# Test the function eosRhoPU.
#
testeosrhopu: testeosrhopu.o $(objects)
	cc -o testeosrhopu testeosrhopu.o $(objects) -lm

#
# Test the function eosRhoPU.
#
testeosurhop: testeosurhop.o $(objects)
	cc -o testeosurhop testeosurhop.o $(objects) -lm

#
# Test the function eosLookupU.
#
testeoslookupu: testeoslookupu.o $(objects)
	cc -o testeoslookupu testeoslookupu.o $(objects) -lm

#
# Test the function tillPressure.
#
testtillpressure: testtillpressure.o $(objects)
	cc -o testtillpressure testtillpressure.o $(objects) -lm

#
# Test the function tillPressureSound.
#
testtillsound: testtillsound.o $(objects)
	cc -o testtillsound testtillsound.o $(objects) -lm

#
# Test the function tillRhoPU().
#
testtillrhopu: testtillrhopu.o $(objects)
	cc -o testtillrhopu testtillrhopu.o $(objects) -lm

#
# Test the ideal gas EOS.
#
testtillidealgas: testtillidealgas.o $(objects)
	cc -o testtillidealgas testtillidealgas.o $(objects) -lm

#
# Generate an ASCII file containing the cold curve u_cold(rho).
#
tillmakecoldcurve: tillmakecoldcurve.o $(objects)
	cc -o tillmakecoldcurve tillmakecoldcurve.o $(objects)  -lm

#
# Test tillColdULookup().
#
testtillcoldu: testtillcoldu.o $(objects)
	cc -o testtillcoldu testtillcoldu.o $(objects)  -lm

#
# Print the Tillotson EOS parameters for a material.
#
tillprintmat: tillprintmat.o $(objects)
	cc -o tillprintmat tillprintmat.o $(objects)  -lm

#
# Test the Woolfson (2007) density correction.
#
testwoolfson: testwoolfson.o woolfson.o $(objects)
	cc -o testwoolfson testwoolfson.o woolfson.o $(objects) -lm

#
# Calculate the Woolfson (2007) coefficients for two given material.
#
calc_fij: calc_fij.o woolfson.o $(objects)
	cc -o calc_fij calc_fij.o woolfson.o $(objects) -lm -Wall

#
# Set a minimum sound speed.
#
soundspeed_cutoff: soundspeed_cutoff.o $(objects)
	cc -o soundspeed_cutoff soundspeed_cutoff.o $(objects) -lm

#
# Calculate P(rho, T=const) for different temperatures.
#
tillpressrhotemp: tillpressrhotemp.o $(objects)
	cc -o tillpressrhotemp tillpressrhotemp.o $(objects) -lm

#
# Check if P(rho, T=const) is monotonic.
#
tillpressrhotemp_monotonic: tillpressrhotemp_monotonic.o $(objects)
	cc -o tillpressrhotemp_monotonic tillpressrhotemp_monotonic.o $(objects) -lm
clean:
	rm $(objects)

cleanall:
	rm $(exe) $(objects) *.o
