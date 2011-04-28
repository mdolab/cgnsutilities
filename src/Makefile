# TECIO LIB
include ../Makefile.in
TECIO = tecio.a

default: all

time_combine:
	$(FC) -c $(FFLAGS) -I$(CGNS_INCLUDE) time_combine.f90
	$(FC) time_combine.o $(TECIO) -o time_combine $(CGNS_LIB) -lstdc++
ts_combine:
	$(FC) -c $(FFLAGS) -I$(CGNS_INCLUDE) ts_combine.f90
	$(FC) ts_combine.o $(TECIO) -o ts_combine $(CGNS_LIB) -lstdc++

cgns_scale:
	$(FC) -c $(FFLAGS) -I$(CGNS_INCLUDE) cgns_scale.f90
	$(FC) cgns_scale.o -o cgns_scale $(CGNS_LIB)

cgns_split:
	$(FC) -c $(FFLAGS) -I$(CGNS_INCLUDE) cgns_split.f90
	$(FC) cgns_split.o -o cgns_split $(CGNS_LIB)

cgns_create:
	f2py  -c  -I$(CGNS_INCLUDE) $(CGNS_LIB)  -m cgns_create cgns_create.f90


all:	time_combine\
	ts_combine\
	cgns_scale \
	cgns_split\
	cgns_create

clean:
	-rm *.o
	-rm *.mod
	-rm time_combine
	-rm ts_combine
	-rm cgns_split
	-rm cgns_scale

