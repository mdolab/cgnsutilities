# Define the CNGS include directory and linking flags for CGNSlib. We
# can use 3.2.x OR CGNS 3.3+. You must define which version is being
# employed as shown below. We are assuming that HDF5 came from PETSc
# so it is included in ${PETSC_LIB}. Otherwise you will have to
# specify the HDF5 library.

# ----------- CGNS ------------------
# CGNS_VERSION_FLAG=               # for CGNS 3.2.x
CGNS_VERSION_FLAG=-DUSECGNSMODULE  # for CGNS 3.3.x
CGNS_INCLUDE_FLAGS=-I$(CGNS_HOME)/include
CGNS_LINKER_FLAGS=-L$(CGNS_HOME)/lib -lcgns


# Intel Fortran Compiler
# CC = gcc
# CFLAGS = -O2 -fPIC
# FC = ifort
# FFLAGS = -O2 -r8 -g -fPIC ${CGNS_VERSION_FLAG}

# Gfortran compiler
CC = gcc
CFLAGS = -O2 -fPIC
FC = gfortran
FFLAGS= -O2 -fdefault-real-8 -g -fPIC ${CGNS_VERSION_FLAG}

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python3-config # use python-config for python 2
F2PY = f2py