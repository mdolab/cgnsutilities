#!/bin/bash
set -e

# Run tests
testflo -v -n 1 --coverage --coverpkg cgnsutilities

# Check that we can run the command line script
cd $HOME
cgns_utils -h
