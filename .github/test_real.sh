#!/bin/bash
set -e
# There are no tests for cgnsutilities
# we just check that we can run the command line script here
cd $HOME
cgns_utils -h