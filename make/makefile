#!/bin/sh
 CC = gfortran
 FFLAGS= -c -O1 
 FFLAGS2= -O3 
 ODIR=../bin/objs
#FFLAGS= -c -O3 -fno-range-check 
#FFLAGS = -c -O0 -g -Wall -Wextra -pedantic -fbounds-check -fno-range-check -fbacktrace
#FFLAGS2= -O3 -g -Wall -Wextra -pedantic -fbounds-check
_OBJS= stringmod.o types.o math.o globals.o physicalconstants.o reading.o hcl2pot.o potinterpol.o pothxhe.o potcal.o  quadrature.o dmc_drv.o main.o 
_OBJS2= hfdimpot.o

OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))
OBJS2 = $(patsubst %,$(ODIR)/%,$(_OBJS2))
../bin/dyncal :  ${OBJS} ${OBJS2} 
	$(CC) ${OBJS} ${OBJS2} -o   ../bin/qdyn  ${FFLAGS2}

$(ODIR)/%.o : ../src/%.f90
	$(CC) ${FFLAGS} $< -o  $@  -J$(ODIR)

$(ODIR)/%.o : ../src/hfdimpot.f
	$(CC) ${FFLAGS} $< -o  $@  -J$(ODIR)

clean: 
	rm $(ODIR)/*.o 
	rm $(ODIR)/*.mod



