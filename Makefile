goal: makefile.dep
	make sd_v2.x 
# superdyson.x and tr_1rdm.x are deprecated

ACT = sed -e 's/^!\*nq/    /' # Disable quad-math statements

IDS := built on $(shell hostname) at $(shell date)

# ifort 18.0.3.222
# F90 = ifort -I. -xAVX2 -qopenmp -O3 -fp-model precise -assume protect_parens -heap-arrays 32 -warn -warn -assume buffered_io -gen-interfaces nosource -mkl -debug full -debug extended -traceback -cpp -D__BUILD_ID__='"Optimized ifort, $(IDS)"'
# LIBS =

# gfortran
# gfortran 7.5.0 works
# OPT
  F90 = gfortran   -I. -O3 -flto -fprotect-parens -fno-fast-math -march=native -mtune=native  -fstack-arrays -fopenmp -fcx-fortran-rules -fno-realloc-lhs -fbacktrace -g -cpp -D__BUILD_ID__='"Optimized gfortran, $(IDS)"'
# F90 = gfortran   -I. -Og -fprotect-parens -fno-fast-math -march=native -mtune=native  -fstack-arrays -fopenmp -std=gnu -pedantic -Wall -Wno-unused-function -Wno-unused-dummy-argument -fcheck=all -fno-realloc-lhs -ffree-line-length-none -cpp -D__BUILD_ID__='"Debug gfortran, $(IDS)"'
  LIBS = -llapack -lblas -lpthread

MAKEFLAGS = -r

.SUFFIXES:
.SUFFIXES: .f90 .o .x .c

OBJS=accuracy.o constants.o block_determinant.o block_diag.o block_matrices.o \
     dgedi.o dgefa.o gamess_internal.o import_gamess.o lapack.o math.o versions.o \
     os_integral_operators.o printing.o sd_core.o sort_tools.o superdyson.o timer.o tr_1rdm.o


sd_v2.x:	sd_v2.f90 $(OBJS)
	$(F90) -o sd_v2.x sd_v2.f90 $(OBJS) $(LIBS)

superdyson.x:	superdyson_driver.o $(OBJS)
	$(F90) -o superdyson.x superdyson_driver.o $(OBJS) $(LIBS)

tr_1rdm.x:	tr_1rdm_driver.o $(OBJS)
	$(F90) -o tr_1rdm.x tr_1rdm_driver.o $(OBJS) $(LIBS)

dgedi.o:	dgedi.f
	$(F90) -c dgedi.f

dgefa.o:	dgefa.f
	$(F90) -c dgefa.f

.f90.o:
	$(ACT) $(ACT2) $< >preprocess/$<
	$(F90) -c preprocess/$<

clean:
	-rm *.o *.mod superdyson.x tr_1rdm.x makefile.dep

makefile.dep: $(shell echo *.f90)
	./make-depend.sh $^ > $@

#
# Automatically-generated dependencies
#
include makefile.dep
