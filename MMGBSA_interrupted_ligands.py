# Place this file in the output folder before executing
# This script extracts the ligands that have been interrupted during the MMGBSA process

import glob
import sys
# list containing the all files to be processed by MMGBSA
list_of_lig_pdb = glob.glob('lig_pdb/*.pdb')

# list containing the files already processed by the MMGBSA
list_of_lig_conformation = glob.glob('1QCF_cluster2*')

fp_lig_recp = open('complexes_to_remove.txt', 'w')

ligands_to_be_removed = []

for each_rec in list_of_lig_pdb:
    each_rec = each_rec.split('/')[1]
    lig_recp_complex = each_rec.replace('.pdb', '')
    if lig_recp_complex not in list_of_lig_conformation:
        lig_recp_name = lig_recp_complex.split('.')[0]
        if lig_recp_name not in ligands_to_be_removed:
            ligand = lig_recp_name.split('_')[2]
            # list containing the ligands that were interrupted during MMGBSA
            ligands_to_be_removed.append(lig_recp_name)
            fp_lig_recp.write(ligand+'\n')

fp_lig_recp.close()
