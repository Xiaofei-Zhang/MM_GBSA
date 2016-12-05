This repository contains some scripts used to do the MM/PBSA analysis on the docking results.

    get_lig_lib.py: Create the library files for all ligands
    change_rcptr.py: Chage the receptor pdb files (used in docking step) to the amber recognizable pdb files
    mmgbsa_MPI.py: Perform the MM/GBSA analysis
    test.sh: sample Slurm batch job submission script (Using UKY DLX Cluster)
    get_avg_energy_best_5.py: Script to collect all MM/GBSA score

Usage:
(1) Create the ligand list file contains the names of ligands to analysis
(2) Create the receptor list file contains the names of receptors to analysis
(3) Verify the directories that contains the MOL2 files of the ligands
(4) Modify the paths information in get_lig_lib.py as needed
(5) Run get_lig_lib.py to generate the library files for all ligands
(7) Copy Y2P.off and frcmod_y2p to amber_lib directory generated in step (5)
(8) Create the .conf file (one input for mmgbsa_MPI.py) contains the correct path information
(10)Run mmpbsa_MPI.py on 1 ligand and 1 receptor to make sure that all inputs are correct.
     eg: mpiexec -np 1 python mmgbsa_MPI.py lig rcptr test.conf
(11)Prepare the batch job script to run mmpbsa_MPI.py


References:
1. D.A. Case, T.A. Darden, T.E. Cheatham, III, C.L. Simmerling, J. Wang, R.E. Duke, R. Luo, R.C. Walker, W. Zhang, K.M. Merz, B. Roberts, S. Hayik, A. Roitberg, G. Seabra, J. Swails, A.W. Götz, I. Kolossváry, K.F. Wong, F. Paesani, J. Vanicek, R.M. Wolf, J. Liu, X. Wu, S.R. Brozell, T. Steinbrecher, H. Gohlke, Q. Cai, X. Ye, J. Wang, M.-J. Hsieh, G. Cui, D.R. Roe, D.H. Mathews, M.G. Seetin, R. Salomon-Ferrer, C. Sagui, V. Babin, T. Luchko, S. Gusarov, A. Kovalenko, and P.A. Kollman (2012), AMBER 12, University of California, San Francisco.

2. MMPBSA.py: An Efficient Program for End-State Free Energy Calculations J. Chem. Theory Comput. In Journal of Chemical Theory and Computation, Vol. 8, No. 9. (25 July 2012), pp. 3314-3321, doi:10.1021/ct300418h by Bill R. Miller, T. Dwight McGee, Jason M. Swails, et al.

3. Parameters for TYR-PO3 : Y2P; N.Homeyer, A.H.C.Horn, H.Lanig, H.Sticht, J. Mol. Model., in press

4. VinaMPI: facilitating multiple receptor high-throughput virtual docking on high-performance computers. Sally R. Ellingson, Jeremy C. Smith, Jerome Baudry J Comput Chem. 2013 Sep 30; 34(25): 2212–2221. doi: 10.1002/jcc.23367
