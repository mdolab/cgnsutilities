# TECIO LIB
include ../Makefile.in
TECIO = tecio.a

default: all

cgns_surf_convert:
	$(FC) -c $(FFLAGS) -I$(CGNS_INCLUDE) cgns_surf_convert.f90
	$(FC) cgns_surf_convert.o $(TECIO) -o cgns_surf_convert $(CGNS_LIB) -lstdc++

time_combine:
	$(FC) -c $(FFLAGS) -I$(CGNS_INCLUDE) time_combine.f90
	$(FC) time_combine.o $(TECIO) -o time_combine $(CGNS_LIB) -lstdc++
ts_combine:
	$(FC) -c $(FFLAGS) -I$(CGNS_INCLUDE) ts_combine.f90
	$(FC) ts_combine.o $(TECIO) -o ts_combine $(CGNS_LIB) -lstdc++

cgns_scale:
	$(FC) -c $(FFLAGS) -I$(CGNS_INCLUDE) cgns_scale.f90
	$(FC) cgns_scale.o -o cgns_scale $(CGNS_LIB)

cgns_flip:
	$(FC) -c $(FFLAGS) -I$(CGNS_INCLUDE) cgns_flip.f90
	$(FC) cgns_flip.o -o cgns_flip $(CGNS_LIB)

cgns_mirror:
	$(FC) -c $(FFLAGS) -I$(CGNS_INCLUDE) cgns_mirror.f90
	$(FC) cgns_mirror.o -o cgns_mirror $(CGNS_LIB)

cgns_divide:
	$(FC) -c $(FFLAGS) -I$(CGNS_INCLUDE) cgns_divide.f90
	$(FC) cgns_divide.o -o cgns_divide $(CGNS_LIB)

cgns_split:
	$(FC) -c $(FFLAGS) -I$(CGNS_INCLUDE) cgns_split.f90
	$(FC) cgns_split.o -o cgns_split $(CGNS_LIB)
cgns_coarsen:
	$(FC) -c $(FFLAGS) -I$(CGNS_INCLUDE) cgns_coarsen.f90
	$(FC) cgns_coarsen.o -o cgns_coarsen $(CGNS_LIB)

cgns_refine:
	$(FC) -c $(FFLAGS) -I$(CGNS_INCLUDE) cgns_refine.f90
	$(FC) cgns_refine.o -o cgns_refine $(CGNS_LIB)

cgns_create:
	f2py  -c --f90flags='$(FFLAGS)' -I$(CGNS_INCLUDE) $(CGNS_LIB)  -m cgns_create cgns_create.f90


surface_extractor:
	$(FC) -c $(FFLAGS) -I$(CGNS_INCLUDE) surface_extractor.f90
	$(FC) surface_extractor.o -o surface_extractor $(CGNS_LIB)


all:	cgns_scale \
	cgns_flip \
	cgns_mirror \
	cgns_split\
	cgns_create\
	cgns_coarsen\
	cgns_refine\
	cgns_divide\
	time_combine\
	ts_combine\
	cgns_surf_convert\
	surface_extractor
clean:
	-rm *.o
	-rm *.mod
	-rm time_combine
	-rm ts_combine
	-rm cgns_split
	-rm cgns_scale
	-rm cgns_flip
	-rm cgns_flip
	-rm cgns_coarsen
	-rm cgns_divide
	-rm cgns_refine
	-rm cgns_mirror
