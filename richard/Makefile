#F90=f95
F90=gfortran
CC=gcc
F77=gfortran

# standard fast flags
FFLAGS= -O3 #-r8

RUNDIR=/local/spg/RUNS/.
INITDIR=/local/spg/RUNS/INITIALS/.

#PG_FLAGS=-L/star/lib -lpgplot -lgwm -L/usr/X11R6/lib -lX11 -lemsf -lems -lcnf -lgcc -lc -lstarmem
#PG_FLAGS=-L/usr/X11R6/lib -lX11 -Wl,-framework -Wl,Foundation -L/sw/lib -lpng -lz -laquaterm -L/sw/lib/pgplot -lpgplot
##PG_FLAGS=-L/usr/lib64 -lpng -lpgplot
PG_FLAGS=-L/usr/lib -L/opt/local/lib -lpgplot #-lpng

.SUFFIXES: .f90 .F90 .f
.f90.o:
	$(F90) $(FFLAGS) -c $<
.F90.o:
	$(F90) $(FFLAGS) $(CFLAGS) -c $<
.f.o:
	$(F77) -c $<

multiple_system: analysis_modules.o multiple_system.o get_centre_of_mass.o get_centre_of_velocity.o binary_orbit.o
	$(F90) $(FFLAGS) -o multiple_system analysis_modules.o multiple_system.o get_centre_of_mass.o get_centre_of_velocity.o binary_orbit.o $(PG_FLAGS)

