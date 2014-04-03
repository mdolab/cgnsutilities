include ../config.mk

PYTHON_OBJECTS = fortranobject.o \
		 libcgns_utilsmodule.o

default:

# Compile the source files:
	$(FC) -c $(FFLAGS) $(CGNS_INCLUDE_FLAGS) cgns_utilities.F90

# Generate Python inlude directory
	 $(eval PYTHON_INCLUDES = $(shell $(PYTHON-CONFIG) --includes))
	 @echo "#------------------------------------------------------#"
	 @echo Python Inclue Flags $(PYTHON_INCLUDES)
	 @echo "#------------------------------------------------------#"

# Generate Numpy inlude directory
	$(eval NUMPY_INCLUDES = $(shell $(PYTHON) -c 'import numpy; print numpy.get_include()'))
	@echo "#------------------------------------------------------#"
	@echo Numpy Include Directory: $(NUMPY_INCLUDES)
	@echo "#------------------------------------------------------#"

# Generate f2py root directory
	$(eval F2PY_ROOT = $(shell $(PYTHON) get_f2py.py))
	@echo "#------------------------------------------------------#"
	@echo f2py root directory: $(F2PY_ROOT)
	@echo "#------------------------------------------------------#"

# Generate pyf file automatically
	f2py -m  --overwrite-signature libcgns_utils -h cgns_utils.pyf cgns_utilities.F90

# Generate wrapper file:
	f2py cgns_utils.pyf

# Compile c wrapper
	$(CC) $(CFLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-I$(F2PY_ROOT)/src -c libcgns_utilsmodule.c

# Compile fortranobject needed by all f2py modules
	$(CC) $(CFLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-c $(F2PY_ROOT)/src/fortranobject.c -o fortranobject.o

# Final Link:
	$(FC) -shared $(PYTHON_OBJECTS) $(FFLAGS) $(CGNS_LINKER_FLAGS) cgns_utilities.o -o libcgns_utils.so
	python importTest.py
	mv libcgns_utils.so ../bin
