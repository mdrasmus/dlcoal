#!/bin/sh
# commands for cleaning up example runs
#
# These commands and meant to be copied and pasted into the command line.
# 


#=============================================================================
# clean up reconciliation files

# remove all MPR reconciliations
find data/flies -name '*mpr*' | xargs rm

# remove all DLCoal reconciliations
find data/flies -name '*dlcoal*' | xargs rm

# clean up simulations
rm -rf data/flies

# clean up relations
rm -rf data/flies-rel


