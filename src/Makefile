include ../config.mk

PYTHON_OBJECTS = fortranobject.o \
		 kd_tree.o \
		 libcgns_utilsmodule.o \
		 libcgns_utils-f2pywrappers2.o


default:

# Compile the source files:
	$(FC) -c $(FFLAGS) $(CGNS_INCLUDE_FLAGS) kd_tree.f90
	$(FC) -c $(FFLAGS) $(CGNS_INCLUDE_FLAGS) -I./ cgns_utilities.F90

# Generate Python include directory
	 $(eval PYTHON_INCLUDES = $(shell $(PYTHON-CONFIG) --includes))
	 @echo "#------------------------------------------------------#"
	 @echo Python Inclue Flags $(PYTHON_INCLUDES)
	 @echo "#------------------------------------------------------#"

# Generate Numpy include directory
	$(eval NUMPY_INCLUDES = $(shell $(PYTHON) -c 'from __future__ import print_function; import numpy; print(numpy.get_include())'))
	@echo "#------------------------------------------------------#"
	@echo Numpy Include Directory: $(NUMPY_INCLUDES)
	@echo "#------------------------------------------------------#"

# Generate f2py root directory
	$(eval F2PY_ROOT = $(shell $(PYTHON) get_f2py.py))
	@echo "#------------------------------------------------------#"
	@echo f2py root directory: $(F2PY_ROOT)
	@echo "#------------------------------------------------------#"

# Generate pyf file automatically
#       f2py -m  --overwrite-signature libcgns_utils -h cgns_utils.pyf cgns_utilities.F90

# Generate wrapper file:
	$(F2PY) cgns_utils.pyf

# Compile c wrapper
	$(CC) $(CFLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-I$(F2PY_ROOT)/src -c libcgns_utilsmodule.c

# Compile the module warapper
	$(FC) $(FFLAGS) -I./ -c libcgns_utils-f2pywrappers2.f90

# Compile fortranobject needed by all f2py modules
	$(CC) $(CFLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-c $(F2PY_ROOT)/src/fortranobject.c -o fortranobject.o

# Final Link:
	$(FC) -shared $(PYTHON_OBJECTS) cgns_utilities.o tecio.a $(CGNS_LINKER_FLAGS) -lstdc++ -o libcgns_utils.so
	$(PYTHON) importTest.py
	mv libcgns_utils.so ../bin
