#!/bin/bash

# This scripts removes the protein-ligand complexes that have been interrupted during the MMGBSA process
# Place this script in the output folder before executing

for f in $(cat complexes_to_remove.txt);
do
rm -r *"$f"*;
done
