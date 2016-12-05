This repository contains some scripts used to do the MM/PBSA analysis on the docking results.

Usage:
(1) Create the ligand list file contains the names of ligands to analysis
(2) Create the receptor list file contains the names of receptors to analysis
(3) Verify the directories that contains the MOL2 files of the ligands
(4) Modify the paths in get_lig_lib.py as needed
(5) Run get_lig_lib.py to generate the library files for all ligands
(7) Copy Y2P.off and frcmod_y2p to amber_lib directory generated in step (5)
(8) Verify the directories that contain the following kinds of files:
	PDBQT files: the docking results generated from VinaMPI (the file name should be in 'receptor_ligand.pdbqt' format)
	Ligand PDB files: the pdb files of ligands in the ligand list file
	Receptor PDB files: the pdb files of receptors in the receptor list file
	Ligand MOL2 files: the directory contains Amber style MOL2 files generated in step (5)
	Amber inputs files: .in files used in Amber molecular dynamics and analysis
(9) Modify the paths in prepare_mmpbsa.py as needed
(10)Run prepare_mmpbsa.py on 1 ligand and 1 receptor to make sure that all paths of files are correct
(11)Prepare the batch job script to run prepare_mmpbsa_MPI.py
