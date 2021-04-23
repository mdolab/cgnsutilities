#!/bin/bash

# This script demonstrates how to run the and use the overwriteFamilies option
CGNS_INFILE="717_wl_L2.cgns"
CGNS_OUTFILE="717_wl_L2_overwriteFamily.cgns"
FAMFILE="overwriteFamily_famFile"
WRITEFILE=true

# Prepare the family file input. Each line has the following format:
# <block number> <faceID> <family>

if [[ $WRITEFILE -eq true ]]; then
cat << 'EOF' > $FAMFILE
1 kLow wing1
EOF
fi

# Run the command
cgns_utils family $CGNS_INFILE $FAMFILE $CGNS_OUTFILE