# cgnsUtilities
[![Build Status](https://dev.azure.com/mdolab/Public/_apis/build/status/mdolab.cgnsutilities?repoName=mdolab%2Fcgnsutilities&branchName=master)](https://dev.azure.com/mdolab/Public/_build/latest?definitionId=16&repoName=mdolab%2Fcgnsutilities&branchName=master)

This repository contains a single program called `cgns_utils` that
provides many useful functions for working with cgns grids. 

## Installation Instructions

First, copy the default configuration file to the required name. For example, to use gfortran: 
  
    cp config/defaults/config.LINUX_GFORTRAN.mk config/config.mk

Open and edit the config file.
You many have to adjust the `CGNS_INCLUDE_FLAGS` and `CGNS_LINKER_FLAGS` to match your installation of the CGNS library.
You may also specify the c-compiler with the `CC` variable and the flags for the c-compiler with `CFLAGS`.
The C-compiler is only used for the compiler `f2py` wrapper.
The Fortran compiler may be specified with the `FC` variable and the corresponding flags with the `FFLAGS` variable.
It has been tested with both Intel and GNU Fortran compilers. 

To compile, simply type

    make 

If all goes well, you will see

    Testing if module libcgns_utils can be imported...
    Module libcgns_utils was successfully imported.

which indicates a successful compilation. 

To install, type

    pip install .

or optionally with the `--user` flag if you are not using a virtual environment.
A console script called `cgns_utils` is provided, which should be installed automatically and available without modifying your `$PATH`.

## Usage

All of the `cgnsUtilties` functionality is accessible through the `cgns_utils` command.
To see a list of the available sub-commands, type

    cgns_utils -h

To get further help for a sub-command, type for example:

    cgns_utils scale -h

## License

Copyright 2020 MDO Lab. See the LICENSE file for details
