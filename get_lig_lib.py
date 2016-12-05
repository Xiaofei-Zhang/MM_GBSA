# get_lig_lib.py
# Build lig.frcmod and lig.lib of the ligands to generate the cmplx.prmtop and 
# cmplx.inpcrd
# Usage:
#	python get_lig_lib $LIGLIST $LIGMOL2DIR
# 
# $LIGLIST is the file contains the names of the ligands
# $LIGMOL2DIR is the directory contains the ligand.mol2 files
#
# Author: Xiaofei Zhang
# Date: July 13 2016

import os
import sys

if len(sys.argv) != 3:
        print("Usage: $LIGLIST $LIGMOL2DIR")
	sys.exit(-1)

# Modify the path of amber/bin and tleap as needed		
amber_bin = '/share/cluster/RHEL6.2/x86_64/apps/amber/openmpi/12p12/bin'
tleap = '/share/cluster/RHEL6.2/x86_64/apps/amber/openmpi/12p12/AmberTools/src/leap/tleap'

if not os.path.isdir('./amber_lib'):
	os.system('mkdir amber_lib')
if not os.path.isdir('./amber_mol2'):
        os.system('mkdir amber_mol2')
if not os.path.isdir('./tleap_inputs'):
        os.system('mkdir tleap_inputs')
lig_dir = sys.argv[2]
lig_list = [line.strip() for line in open(sys.argv[1], 'r')]



for lig in lig_list:
	ligbase = lig.split('_')[0]
	lig_file = lig_dir + '/' + lig + '.mol2'
	amber_mol = lig + '_amber.mol2'
	cmd1 = amber_bin + '/antechamber -i ' + lig_file + ' -fi mol2 -o '+ amber_mol +' -fo mol2 -s 2'
	os.system(cmd1)
	
	#renaming residue in amber.mol2 file
	cmd_rn = 'sed -i \'s/<0>/CHE/g; s/'+ligbase+'/CHE/g\' ' + amber_mol
	os.system(cmd_rn)	
	frcmod = lig + '.frcmod'
	cmd2 = amber_bin + '/parmchk -i ' + amber_mol + ' -f mol2 -o ' + frcmod
	os.system(cmd2)
	tleap_fname = lig + '_tleap.in'
	tleap_lib = lig + '.lib'
	with open(tleap_fname, 'w') as tleap_f:
		tleap_f.write('source leaprc.ff12SB\n')
		tleap_f.write('source leaprc.gaff\n')
		tleap_f.write('CHE = loadmol2 '+ amber_mol + '\n')
		tleap_f.write('check CHE\n')
		tleap_f.write('loadamberparams ' + frcmod +'\n')
		tleap_f.write('saveoff CHE ' + tleap_lib +'\n')
		tleap_f.write('saveamberparm CHE ' + lig + '.prmtop ' + lig + '.inpcrd' +'\n')
		tleap_f.write('quit' +'\n')
	cmd3 = tleap + ' -f ' + tleap_fname
	os.system(cmd3)
	os.system('mv ./*.frcmod ./amber_lib/')
	os.system('mv ./*.prmtop ./amber_lib/')
	os.system('mv ./*.inpcrd ./amber_lib/')
	os.system('mv ./*.lib ./amber_lib/')
	os.system('mv ./*.mol2 ./amber_mol2/')
	os.system('mv ./*.in ./tleap_inputs/')


