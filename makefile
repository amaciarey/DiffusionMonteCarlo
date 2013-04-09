f90=ifort -fast
#f90=gfortran -Wall -fno-range-check -C

dmc: dmc.o\
 	dmc_mod.mod\
 	dmc_mod.o\
 	random_mod.mod\
 	random_mod.o\
 	sample_mod.o\
 	sample_mod.mod\
 	system_mod.mod\
 	system_mod.o\
 	r8_gamma.o\
 	interpolate.o
	$(f90) -o dmc *.o

dmc.o: dmc.f90\
 	global_mod.mod\
 	random_mod.mod\
 	sample_mod.mod\
 	dmc_mod.mod
	$(f90) -c dmc.f90

dmc_mod.o dmc_mod.mod: dmc_mod.f90\
 	global_mod.mod\
 	random_mod.mod\
 	system_mod.mod\
 	sample_mod.mod\
 	pbc_mod.mod\
 	interpolate.f90
	$(f90) -c dmc_mod.f90

sample_mod.o sample_mod.mod: sample_mod.f90\
 	global_mod.mod\
 	random_mod.mod\
 	system_mod.mod\
 	pbc_mod.mod
	$(f90) -c sample_mod.f90

system_mod.o system_mod.mod: system_mod.f90\
 	global_mod.mod\
 	bessel_mod.mod
	$(f90) -c system_mod.f90

pbc_mod.o pbc_mod.mod: pbc_mod.f90\
 	global_mod.mod
	$(f90) -c pbc_mod.f90

bessel_mod.o bessel_mod.mod: bessel_mod.f90
	$(f90) -c bessel_mod.f90

global_mod.o global_mod.mod: global_mod.f90
	$(f90) -c global_mod.f90

random_mod.o random_mod.mod: random_mod.f90
	$(f90) -c random_mod.f90

r8_gamma.o: r8_gamma.f90
	$(f90) -c r8_gamma.f90

interpolate.o: interpolate.f90
	$(f90) -c interpolate.f90

clean:
	rm -f *.o *.mod 
