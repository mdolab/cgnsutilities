#!/bin/bash
set -e

# This script demonstrates how to run and use the overwritebc option to add array bc data
CGNS_INFILE="717_wl_L2.cgns"
CGNS_OUTFILE="717_wl_L2_overwriteBCs_array.cgns"
BC_FILE="hotwall_boco.info"


# Run the command
cgns_utils overwriteBC $CGNS_INFILE $BC_FILE $CGNS_OUTFILE
