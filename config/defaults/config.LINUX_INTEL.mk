# ------- Define CGNS Inlcude and linker flags -------------------------
# Define the CGNS include directory and linking flags for the CGNS library.
# We are assuming that HDF5 came from PETSc so it is included in ${PETSC_LIB}.
# Otherwise you will have to specify the HDF5 library.
CGNS_INCLUDE_FLAGS=-I$(CGNS_HOME)/include
CGNS_LINKER_FLAGS=-L$(CGNS_HOME)/lib -lcgns

# Intel compilers
ICC_EXISTS := $(shell command -v icc)
ifdef ICC_EXISTS
  # icc only exists on older Intel versions
  # Assume that we want to use the old compilers
  FC = ifort
  CC = icc
else
  # Use the new compilers
  FC = ifx
  CC = icx
endif

# Compiler flags
CFLAGS = -O2 -fPIC
FFLAGS = -O2 -r8 -g -fPIC -stand f08

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python3-config # use python-config for python 2
F2PY = f2py
