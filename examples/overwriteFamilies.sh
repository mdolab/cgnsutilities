#!/bin/bash

# This script demonstrates how to run and use the family option
CGNS_INFILE="717_wl_L2.cgns"
CGNS_OUTFILE="717_wl_L2_overwriteFamilies.cgns"
FAMFILE="family_famFile"
WRITEFILE=true

# Prepare the family file input. Each line has the following format:
# <block number> <faceID> <family>
# If using a volume mesh, then <block number> is the domain number from tecplot, not the zone number

if [[ $WRITEFILE -eq true ]]; then
cat << 'EOF' > $FAMFILE
1 kLow wing1
EOF
fi

# Run the command
cgns_utils overwriteFamilies $CGNS_INFILE $FAMFILE $CGNS_OUTFILE
