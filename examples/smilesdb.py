#######################
# Script smilesdb.py  #
# Release: All        #
#######################

# Test database script operation - return a SMILES string on input regno 12345

import sys

regno = '12345'
smiles_string = 'Cc1ncc(COP(O)(O)=O)c(C=O)c1O'

number_of_args = len(sys.argv)

if number_of_args == 2:
    if sys.argv[1] == regno:
        print smiles_string
        
