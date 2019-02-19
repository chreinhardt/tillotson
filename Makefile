objects = tillotson.o tillinitlookup.o tillsplint.o interpol/brent.o 

exe = table pressureoldnew lookup lookup_cold testu1 testspline testsplint testnewsplint testsplint2 testsplinerho testsplintrho testsplinev testsplintv testcubicintrho testlookupucold testudrho testudv testgrid testpolyv printderiv printpress pressneg testisintable testisbelowcoldcurve testrhomin testoutofbounds testsolvebc calcisentrope testrhoptemp tillcalcpressure testdirectintegration testpoverrho2 testeospressure testtillpressure testtillrhopu tillpressrhotemp

#defs = -DTILL_PRESS_NP -DTILL_OUTPUT_ALL_WARNINGS -DTILL_PRESS_MELOSH
defs = -DTILL_PRESS_NP -DTILL_OUTPUT_ALL_WARNINGS -DTILL_VERBOSE

# GNU Science library (uncomment if not needed)
GSL_LIB = -lgsl -lgslcblas

CFLAGS ?= -O3 $(defs) -Wall

LIBS ?= -lm $(GSL_LIB)

default:
	@echo "Please specify which tool you want to make."

all:
	default

# Debug the function tillLookupU() by comparing the results of the lookup to a
# direct calculation using tillCalcU() which simply solves the ODE using a RK4
# integrator.
table: table.o $(objects)
	cc -o table table.o $(objects) $(LIBS)

# Compare the results of tillPressureSoundold() that contains the old code we used in
# Gasoline and tillPressureSound()
pressureoldnew: pressureoldnew.o $(objects)
	cc -o pressureoldnew pressureoldnew.o $(objects) $(LIBS)

#
# Generate a lookup table for a given material.
#
lookup: lookup.o $(objects)
	cc -o lookup lookup.o $(objects) $(LIBS)

#
# Make a lookup table for the cold curve.
#
lookup_cold: lookup_cold.o $(objects)
	cc -o lookup_cold lookup_cold.o $(objects) $(LIBS)

testu1: testu1.o $(objects)
	cc -o testu1 testu1.o $(objects) $(LIBS)

testspline: testspline.o $(objects)
	cc -o testspline testspline.o $(objects) $(LIBS)

#
# Code for debugging the interpolation function tillCubicInt().
#
testsplint: testsplint.o $(objects)
	cc -o testsplint testsplint.o $(objects) $(LIBS)

#
# Code for debugging the new logarithmic lookup table.
#
testsplintlogrho: testsplintlogrho.o $(objects)
	cc -o testsplintlogrho testsplintlogrho.o $(objects) $(LIBS)

#
# Code for debugging tillLookupU using the new logarithmic lookup table.
#
testlookupulogrho: testlookupulogrho.o $(objects)
	cc -o testlookupulogrho testlookupulogrho.o $(objects) $(LIBS)

#
# Code for debugging tillCalcU using log(rho) as integration variable.
#
testcalcu: testcalcu.o $(objects)
	cc -o testcalcu testcalcu.o $(objects) $(LIBS)

#
# Pretty much the same but it compares the old with the new interpolator.
#
testnewsplint: testnewsplint.o $(objects)
	cc -o testnewsplint testnewsplint.o $(objects) $(LIBS)

testsplint2: testsplint2.o $(objects)
	cc -o testsplint2 testsplint2.o $(objects) $(LIBS)

testsplinerho: testsplinerho.o $(objects)
	cc -o testsplinerho testsplinerho.o $(objects) $(LIBS)

testsplintrho: testsplintrho.o $(objects)
	cc -o testsplintrho testsplintrho.o $(objects) $(LIBS)

testsplinev: testsplinev.o $(objects)
	cc -o testsplinev testsplinev.o $(objects) $(LIBS)

testsplintv: testsplintv.o $(objects)
	cc -o testsplintv testsplintv.o $(objects) $(LIBS)

testcubicintrho: testcubicintrho.o $(objects)
	cc -o testcubicintrho testcubicintrho.o $(objects) $(LIBS)

testlookupucold: testlookupucold.o $(objects)
	cc -o testlookupucold testlookupucold.o $(objects) $(LIBS)

testudrho: testudrho.o $(objects)
	cc -o testudrho testudrho.o $(objects) $(LIBS)

testudv: testudv.o $(objects)
	cc -o testudv testudv.o $(objects) $(LIBS)

testgrid: testgrid.o $(objects)
	cc -o testgrid testgrid.o $(objects) $(LIBS)

testpolyv: testpolyv.o $(objects)
	cc -o testpolyv testpolyv.o $(objects) $(LIBS)

printderiv: printderiv.o $(objects)
	cc -o printderiv printderiv.o $(objects) $(LIBS)

printpress: printpress.o $(objects)
	cc -o printpress printpress.o $(objects) $(LIBS)

pressneg: pressneg.o $(objects)
	cc -o pressneg pressneg.o $(objects) $(LIBS)

testisintable: testisintable.o $(objects)
	cc -o testisintable testisintable.o $(objects) $(LIBS)

testisbelowcoldcurve: testisbelowcoldcurve.o $(objects)
	cc -o testisbelowcoldcurve testisbelowcoldcurve.o $(objects) $(LIBS)

testrhomin: testrhomin.o $(objects)
	cc -o testrhomin testrhomin.o $(objects) $(LIBS)

testoutofbounds: testoutofbounds.o $(objects)
	cc -o testoutofbounds testoutofbounds.o $(objects) $(LIBS)

testsolvebc: testsolvebc.o $(objects)
	cc -o testsolvebc testsolvebc.o $(objects) $(LIBS)

calcisentrope: calcisentrope.o $(objects)
	cc -o calcisentrope calcisentrope.o $(objects) $(LIBS)

testrhoptemp: testrhoptemp.o $(objects)
	cc -o testrhoptemp testrhoptemp.o $(objects) $(LIBS)

tillcalcpressure: tillcalcpressure.o $(objects)
	cc -o tillcalcpressure tillcalcpressure.o $(objects) $(LIBS)

#
# Calculate the bulk modulus for different isentropes.
#
tillcalcbulkmodulus: tillcalcbulkmodulus.o $(objects)
	cc -o tillcalcbulkmodulus tillcalcbulkmodulus.o $(objects) $(LIBS)

# Code to test the direct integration tillCalcU().
testdirectintegration: testdirectintegration.o $(objects)
	cc -o testdirectintegration testdirectintegration.o $(objects) $(LIBS)

#
# Print PoverRho2 for small rho to see, if the expression diverges.
#
testpoverrho2: testpoverrho2.o $(objects)
	cc -o testpoverrho2 testpoverrho2.o $(objects) $(LIBS)

#
# Test the function eosPressureSound.
#
testeospressure: testeospressure.o $(objects)
	cc -o testeospressure testeospressure.o $(objects) $(LIBS)

#
# Test the function eosPressureSoundRhoT.
#
testeospressurerhot: testeospressurerhot.o $(objects)
	cc -o testeospressurerhot testeospressurerhot.o $(objects) $(LIBS)

#
# Test the function eosdPdrho.
#
testeosdpdrho: testeosdpdrho.o $(objects)
	cc -o testeosdpdrho testeosdpdrho.o $(objects) $(LIBS)

#
# Test the function eosdPdu.
#
testeosdpdu: testeosdpdu.o $(objects)
	cc -o testeosdpdu testeosdpdu.o $(objects) $(LIBS)

#
# Test the function eosRhoPU.
#
testeosrhopu: testeosrhopu.o $(objects)
	cc -o testeosrhopu testeosrhopu.o $(objects) $(LIBS)

#
# Test the function eosRhoPU.
#
testeosurhop: testeosurhop.o $(objects)
	cc -o testeosurhop testeosurhop.o $(objects) $(LIBS)

#
# Test the function eosLookupU.
#
testeoslookupu: testeoslookupu.o $(objects)
	cc -o testeoslookupu testeoslookupu.o $(objects) $(LIBS)

#
# Test the function tillPressure.
#
testtillpressure: testtillpressure.o $(objects)
	cc -o testtillpressure testtillpressure.o $(objects) $(LIBS)

#
# Test the function tillPressureSound.
#
testtillsound: testtillsound.o $(objects)
	cc -o testtillsound testtillsound.o $(objects) $(LIBS)

#
# Test the function tillRhoPU().
#
testtillrhopu: testtillrhopu.o $(objects)
	cc -o testtillrhopu testtillrhopu.o $(objects) $(LIBS)

#
# Test the ideal gas EOS.
#
testtillidealgas: testtillidealgas.o $(objects)
	cc -o testtillidealgas testtillidealgas.o $(objects) $(LIBS)

#
# Generate an ASCII file containing the cold curve u_cold(rho).
#
tillmakecoldcurve: tillmakecoldcurve.o $(objects)
	cc -o tillmakecoldcurve tillmakecoldcurve.o $(objects)  $(LIBS)

#
# Test tillColdULookup().
#
testtillcoldu: testtillcoldu.o $(objects)
	cc -o testtillcoldu testtillcoldu.o $(objects)  $(LIBS)

#
# Calculate the cold energy for a given density.
#
till_calc_u_cold: till_calc_u_cold.o $(objects)
	cc -o till_calc_u_cold till_calc_u_cold.o $(objects)  $(LIBS)

#
# Print the Tillotson EOS parameters for a material.
#
tillprintmat: tillprintmat.o $(objects)
	cc -o tillprintmat tillprintmat.o $(objects)  $(LIBS)

#
# Test the Woolfson (2007) density correction.
#
testwoolfson: testwoolfson.o woolfson.o $(objects)
	cc -o testwoolfson testwoolfson.o woolfson.o $(objects) $(LIBS)

#
# Calculate the Woolfson (2007) coefficients for two given material.
#
calc_fij: calc_fij.o woolfson.o $(objects)
	cc -o calc_fij calc_fij.o woolfson.o $(objects) $(LIBS) -Wall

#
# Set a minimum sound speed.
#
soundspeed_cutoff: soundspeed_cutoff.o $(objects)
	cc -o soundspeed_cutoff soundspeed_cutoff.o $(objects) $(LIBS)

#
# Calculate P(rho, T=const) for different temperatures.
#
tillpressrhotemp: tillpressrhotemp.o $(objects)
	cc -o tillpressrhotemp tillpressrhotemp.o $(objects) $(LIBS)

#
# Check if P(rho, T=const) is monotonic.
#
tillpressrhotemp_monotonic: tillpressrhotemp_monotonic.o $(objects)
	cc -o tillpressrhotemp_monotonic tillpressrhotemp_monotonic.o $(objects) $(LIBS)

clean:
	rm $(objects)

cleanall:
	rm $(exe) $(objects) *.o
