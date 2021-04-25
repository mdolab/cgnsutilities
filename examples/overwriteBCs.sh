#!/bin/bash

# This script demonstrates how to run and use the overwritebc option
CGNS_INFILE="717_wl_L2.cgns"
CGNS_OUTFILE="717_wl_L2_overwriteBCs.cgns"
BC_FILE="overwriteBCs_bcFile"
WRITE_FILE=true

# Prepare the BC input file we want to use. Each line has the following format:
# <block number> <faceID> <BCType> <family> [dataset]

if [[ $WRITE_FILE -eq true ]]; then
cat << 'EOF' > $BC_FILE
1 kLow bcwallinviscid wall_inviscid
2 kLow bcoutflowsubsonic wing_le_inlet BCDataSet_1 BCInFlowSubsonic Dirichlet PressureStagnation 1234.0 TemperatureStagnation 4556.0
EOF
fi

# Run the command
cgns_utils overwritebc $CGNS_INFILE $BC_FILE $CGNS_OUTFILE