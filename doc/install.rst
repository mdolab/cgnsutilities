Installation Instructions
=========================

First, copy the default configuration file to the required name. For example, to use gfortran:

.. prompt:: bash

    cp config/defaults/config.LINUX_GFORTRAN.mk config/config.mk

Open and edit the config file.
You many have to adjust the ``CGNS_INCLUDE_FLAGS`` and ``CGNS_LINKER_FLAGS`` to match your installation of the CGNS library.
You may also specify the C compiler with the ``CC`` variable and the flags for the C compiler with ``CFLAGS``.
The C-compiler is only used for the compiler ``f2py`` wrapper.
The Fortran compiler may be specified with the ``FC`` variable and the corresponding flags with the ``FFLAGS`` variable.
It has been tested with both Intel and GNU Fortran compilers.

To compile, simply type

.. prompt:: bash

    make 

If all goes well, you will see::

    Testing if module libcgns_utils can be imported...
    Module libcgns_utils was successfully imported.

which indicates a successful compilation.

To install, type

.. prompt:: bash

    pip install .

or optionally with the ``--user`` flag if you are not using a virtual environment.
A console script called ``cgns_utils``` is provided, which should be installed automatically and available without modifying your ``$PATH``.