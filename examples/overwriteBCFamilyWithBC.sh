#!/bin/bash
set -e

# This script demonstrates how to run and use the overwriteBCFamilyWithBC option
CGNS_INFILE="717_wl_L2.cgns"
CGNS_OUTFILE="717_wl_L2_overwriteBCFamilyWithBC.cgns"

# Specify the CGNS family we want to update, the new BC type and which blocks to update
FAMILYNAME="Far"
NEWBCTYPE="bcoverset"
BLOCKIDS="2 4"

# Run the command
cgns_utils overwriteBCFamilyWithBC $CGNS_INFILE $FAMILYNAME $NEWBCTYPE $CGNS_OUTFILE --blockIDs $BLOCKIDS
