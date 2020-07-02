# cgnsUtilities

This repository contains a single program called `cgns_utils` that
provides many useful functions for working with cgns grids. 

## Installation Instructions

First, copy the default configuration file to the required name:
  
    cp config.mk.info config.mk

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

This will show the following sub-command options:

    Choose one of the listed operations to perform
    scale               Scale a grid by a constant factor
    flip                Flip a grid about a plane defined by an axis
    coarsen             Coarsen a grid uniformly
    refine              Refine a grid uniformly
    extract             Extract a wall surface from file
    mirror              Mirror a grid about a plane defined by an axis. This
                        doubles the grid size
    split               Face-match a grid. If the grid is already faced
                        matched, this witll have no effect
    connect             Determine the block-to-block connectivity information
                        for a face-matched grid
    divide              Divide all blocks in the grid into 8 sub-blocks
    autobc              Try to determine boundary conditions for blocks. Only
                        suitable for external flow applications.
    family              Overwrite family information
    rebunch             Rebunch offwall spacing (experimental)
    cgns2plot3d         Convert a cgns file to a plain plot3d file


To get further help for a sub-command, type for example:

    cgns_utils scale -h

which will yield the following:

    usage: cgns_utils scale [-h] gridFile scale [outFile]

    positional arguments:
      gridFile    Name of input CGNS file
      scale       scale factor
      outFile     Optional output file


## License

Copyright 2020 MDO Lab. See the LICENSE file for details
