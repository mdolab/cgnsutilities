#!/bin/bash
set -e

# Check that we can run the command line script
cd $HOME
cgns_utils -h

# Run tests
testflo -v -n 1 .
