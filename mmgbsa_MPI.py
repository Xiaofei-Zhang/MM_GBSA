# mmgbsa_MPI.py
# MM/PBSA analysis script in parallel
# Usage:
#       mpiexec -np #processes python mmgbsa_MPI.py $LIGLIST $RCPTRLIST $CONFIG
#
# $LIG_list is the file contains all the ligand name
# $RCPTR_list is the file contains all the receptor name
# $CONFIG is the configuration file that contains the paths of the programs and # files
#
# Author: Xiaofei Zhang
# Date: Nov 2 2016

import os
import sys
import shlex, subprocess
from mpi4py import MPI

# def main
def main():
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    par_dict = {}
    split_job_dict = {}
    docking_results = {}
    mmgbsa_job_dict = {}
    
    # Initialize the parameters
    if (rank == 0):
        if len(sys.argv) != 4:
            print("Usage: mpiexec -np #processes python prepare_mmpbsa.py $LIGLIST $RCPTRLIST $CONFIG")
            MPI.Finalize()
            sys.exit(-1)
        lig_list = [line.strip() for line in open(sys.argv[1], 'r')]
        rcptr_list = [line.strip() for line in open(sys.argv[2], 'r')]
        split_index = 0
        for rcptr in rcptr_list:
            for lig in lig_list:
                split_job_dict[split_index]=[lig, rcptr]
                split_index += 1
        par_dict = parse_config(sys.argv[3])
        if not os.path.isdir(par_dict['out_dir']):
            os.system('mkdir -p ' + par_dict['out_dir'])
        with cd(par_dict['out_dir']):
            if not os.path.isdir('lig_pdb'):
                os.system('mkdir -p lig_pdb')

    par_dict = comm.bcast(par_dict, root=0)
    split_job_dict = comm.bcast(split_job_dict, root=0)
    
    comm.Barrier()
    
    for key in split_job_dict:
        if (int(key) % size == rank):
            docking_results.update(split_model(split_job_dict[key][0], split_job_dict[key][1], par_dict))
    comm.Barrier()
    
    if rank != 0:
        comm.send(docking_results, dest=0, tag=1)
    else:
        for i in range(size - 1):
            temp = comm.recv(source=MPI.ANY_SOURCE, tag=1)
            docking_results.update(temp)
    comm.Barrier()
    
    if rank == 0:
        job_index=0
        for key in docking_results:
            mmgbsa_job_dict[job_index] = key
            job_index += 1
        print 'Total number of mmgbsa jobs is: ' + str(job_index)
    comm.Barrier()
    mmgbsa_job_dict = comm.bcast(mmgbsa_job_dict, root=0)
    comm.Barrier()
    for key in mmgbsa_job_dict:
        if (int(key) % size == rank):
            mmgbsa_pipeline(mmgbsa_job_dict[key], par_dict)
    comm.Barrier()
   # if rank == 0:
   #     print mmgbsa_job_dict
    MPI.Finalize()

# mmgbsa_pipeline: function to do the whole  MM/GBSA analysis for 1 pose of docking results
def mmgbsa_pipeline(job_name, par_dict):
    tokens = job_name.split('.')
    if len(tokens) != 3:
        print("ERROR: Bad format of mmgbsa job name. Make sure no '.' in lig/rcptr name.")
        return 0
    with cd(par_dict['out_dir']):
        basename_dir = tokens[0] + '_' + tokens[1] + '.' +  tokens[2]
        os.system('mkdir -p ' + basename_dir)
        with cd(basename_dir):
            add_hydrogen(tokens[0], tokens[1], tokens[2], par_dict)
            build_cmplx(tokens[0], tokens[1], tokens[2], par_dict)
            minimization(tokens[0], tokens[1], tokens[2], par_dict)
            mmgbsa(tokens[0], tokens[1], tokens[2], par_dict)

# mmgbsa_analysis: function to do the mmgbsa analysis step
def mmgbsa(rcptr_name, lig_name, pose_number, par_dict):
    prmtop = 'cmplx_' + rcptr_name + '_' + lig_name + '.' + pose_number + '.prmtop'
    inpcrd = 'cmplx_' + rcptr_name + '_' + lig_name + '.' + pose_number + '.inpcrd'
    cmd1 = par_dict['cpptraj'] + ' ' + prmtop + ' <' + par_dict['m_in_dir'] + '/'+ par_dict['cmplx_in'] + '> cmplx.out'
    #print cmd1
    os.system(cmd1)

    cmd2 = par_dict['cpptraj'] + ' ' + prmtop + ' <' + par_dict['m_in_dir'] + '/'+ par_dict['rcptr_in'] + '> rcptr.out'
    #print cmd2
    os.system(cmd2)

    cmd3 = par_dict['cpptraj'] + ' ' + prmtop + ' <' + par_dict['m_in_dir'] + '/'+ par_dict['lig_in'] + '> lig.out'
    #print cmd3
    os.system(cmd3)

    cmplx_inp_prmtop = 'cmplx.strip.strip.' + prmtop
    rcptr_inp_prmtop = 'rcptr.strip.strip.strip.' + prmtop
    lig_inp_prmtop = 'lig.strip.strip.strip.' + prmtop
    cmplx_i = 'cmplx.' + prmtop
    rcptr_i = 'rcptr.' + prmtop
    lig_i = 'lig.' + prmtop
    
    change_prmtop(cmplx_inp_prmtop, cmplx_i)
    change_prmtop(rcptr_inp_prmtop, rcptr_i)
    change_prmtop(lig_inp_prmtop, lig_i)
    
    cmd4 = par_dict['mmpbsa'] + ' -O -i ' + par_dict['m_in_dir'] + '/' + par_dict['mmpbsa_in'] + ' -o FINAL_RESULTS_MMPBSA.dat -cp ' + cmplx_i + ' -rp ' + rcptr_i + ' -lp ' + lig_i + ' -y cmplx.dcd'
    #print cmd4
    os.system(cmd4)

# minimization: function to run the minimization
def minimization(rcptr_name, lig_name, pose_number, par_dict):
    prmtop = 'cmplx_' + rcptr_name + '_' + lig_name + '.' + pose_number + '.prmtop'
    inpcrd = 'cmplx_' + rcptr_name + '_' + lig_name + '.' + pose_number + '.inpcrd'
    cmd1 = 'time ' + par_dict['pmemd'] + ' -O -i ' + par_dict['m_in_dir'] + '/'+ par_dict['min_1_in'] + ' -o min1.out -p ' + prmtop + ' -c ' + inpcrd + ' -ref ' + inpcrd +' -r min1.rst -x min1.mdcrd'
    #print cmd1
    os.system(cmd1)
    cmd2 = 'time ' + par_dict['pmemd'] + ' -O -i ' + par_dict['m_in_dir'] + '/'+ par_dict['min_2_in'] + ' -o min2.out -p ' + prmtop + ' -c min1.rst -ref ' + inpcrd +' -r min2.rst -x min2.mdcrd'
    #print cmd2
    os.system(cmd2)

# build_cmplx: build the complex pdb of ligand and receptor, and build the 
# parameter files for the complex
def build_cmplx(rcptr_name, lig_name, pose_number, par_dict):
    rcptr_pdb = par_dict['rcptr_dir'] + rcptr_name + '.pdb'
    lig_h_pdb = rcptr_name + '_' + lig_name + '.' + pose_number + '.h.pdb'
    cmplx_pdb = 'cmplx_' + rcptr_name + '_' + lig_name + '.' + pose_number + '.pdb'
    cmd1 = 'cat '+ rcptr_pdb + ' > ' + cmplx_pdb
    cmd2 = 'cat '+ lig_h_pdb + ' >> ' + cmplx_pdb
    os.system(cmd1)
    os.system(cmd2)
    tleap_fname = 'cmplx_' + rcptr_name + '_' + lig_name + '.' + pose_number + '_tleap.in'
    tleap_lib = lig_name + '.lib'
    frcmod = lig_name + '.frcmod'
    with open(tleap_fname, 'w') as tleap_f:
        tleap_f.write('source leaprc.ff12SB\n')
        tleap_f.write('source leaprc.gaff\n')
        tleap_f.write('loadamberparams ' + par_dict['amber_lib_dir'] +'/frcmod_y2p\n')
        tleap_f.write('loadoff ' + par_dict['amber_lib_dir'] +'/Y2P.off\n')
        tleap_f.write('loadoff ions08.lib\n')
        tleap_f.write('loadamberparams frcmod.ionsjc_tip3p\n')
        tleap_f.write('loadamberparams ' + par_dict['amber_lib_dir'] +'/'+ frcmod + '\n')
        tleap_f.write('loadoff ' + par_dict['amber_lib_dir'] +'/'+ tleap_lib + '\n')
        tleap_f.write('mol = loadpdb '+ cmplx_pdb + '\n')
        tleap_f.write('solvatebox mol TIP3PBOX 20.0\n')
        tleap_f.write('addIonsRand mol K+ 136 Cl- 124 2\n')
        tleap_f.write('saveamberparm mol ' + cmplx_pdb[:-4] + '.prmtop ' + cmplx_pdb[:-4] +'.inpcrd\n')
        tleap_f.write('quit' +'\n')
    cmd3 = par_dict['tleap'] + ' -f ' + tleap_fname + ' > build_cmplx.log'
    os.system(cmd3)


# add_hydrogen: add the hydrogen atoms to the ligand pdb file
def add_hydrogen(rcptr_name, lig_name, pose_number, par_dict):
    lig_pdb_fname = par_dict['out_dir'] +'/lig_pdb/' + rcptr_name + '_' + lig_name + '.' + pose_number + '.pdb' 
    tmp_pdb = rcptr_name + '_' + lig_name + '.' + pose_number + '.tmp'
    lig_h_pdb = rcptr_name + '_' + lig_name + '.' + pose_number + '.h.pdb'
    lig_mol2 = par_dict['lig_dir'] + '/' + lig_name + '_amber.mol2'
    tleap_fname = rcptr_name + '_' + lig_name + '_add_h'+ '.' + pose_number +'.in'
    atom = []
    atom_list = []
    with open(lig_pdb_fname, 'r') as infile:
        for line in infile:
            atom = parse_atom_line(line)
            atom_list.append(atom)
    with open(tmp_pdb, 'w') as tmpfile:
        for myatom in atom_list:
            write_pdb(myatom, tmpfile)
    with open(tleap_fname, 'w') as tleap_f:
        tleap_f.write('CHE = loadmol2 ' + lig_mol2 + '\n')
        tleap_f.write('list\n')
        tleap_f.write('mol1 = loadpdb ' + tmp_pdb + '\n')
        tleap_f.write('savePdb mol1 ' + lig_h_pdb + '\n')
        tleap_f.write('quit' +'\n')
    cmd = par_dict['tleap'] + ' -f ' + tleap_fname + ' > add_hydrogen.log'
    # print cmd
    os.system(cmd)
    os.system('rm ' + tmp_pdb)
    os.system('rm ' + tleap_fname)

# function to split the docking pdbqt file to get the ligands pdb files
def split_model(lig_name, rcptr_name, par_dict):
    basename = rcptr_name + '_' + lig_name
    # using '.' as the seperator of ligand and receptor name in mmgbsa job dictionary
    dictionary_name = rcptr_name + '.' + lig_name
    rcptr_pdb = par_dict['rcptr_dir'] + '/' + rcptr_name + '.pdb'
    pdbqt = par_dict['pdbqt_dir'] + '/' + rcptr_name + '_' + lig_name + '.pdbqt'
    docking_results = {}
    # pdbqt file split
    with cd(par_dict['out_dir']):
        lig_pdb = open(basename + '.test.txt', 'w')
        lig_pdb.close()
        os.system('rm ' + basename + '.test.txt')
        pdbqt_file = open(pdbqt, 'r')
        model_number = int(par_dict['model_number'])
        model_index = 0
        i = 0
        for line in pdbqt_file:
            if(i == 1):
                lig_pdb_fname = 'lig_pdb/' + basename + '.' + str(model_index) + '.pdb'
                energy_line = line.split()
                if len(energy_line) > 3 and model_index < model_number:
                    # using '.' as the seperator of ligand name and model_index
                    docking_results[dictionary_name + '.' + str(model_index)] = float(energy_line[3])
                    #print pdb_fname
                lig_pdb = open(lig_pdb_fname, 'w')
            if(line[0:5]=='MODEL'):
                i=0
                model_index = model_index + 1
                lig_pdb.close()
            if(i > 0):
                if (line[0:4] == 'ATOM'):
                    lig_pdb.write(line)
                    #pdbfile.write('\n')
            i=i+1
        lig_pdb.close()
        pdbqt_file.close()
        return docking_results


# parse_config
# Parse the configuration file. Return the dict contains the parameter information
def parse_config(config_filename):
    # keys are:
    #    pdbqt_dir
    #    lig_dir
    #    rcptr_dir
    #    amber_lib_dir
    #    out_dir
    #    m_in_dir
    #    min_1_in
    #    min_2_in
    #    cmplx_in
    #    rcptr_in
    #    lig_in
    #    mmpbsa_in
    #    tleap
    #    pmemd
    #    cpptraj
    #    mmpbsa
    #    model_number
    par_dict = {
        'pdbqt_dir':'',
        'lig_dir':'',
        'rcptr_dir':'',
        'amber_lib_dir':'',
        'out_dir':os.getcwd() + '/out/',
        'm_in_dir':'',
        'min_1_in':'',
        'min_2_in':'',
        'cmplx_in':'',
        'rcptr_in':'',    'lig_in':'',
        'mmpbsa_in':'',
        'tleap':'',
        'pmemd':'pmemd',
        'cpptraj':'cpptraj',
        'mmpbsa':'MMPBSA.py',
        'model_number':20
    }
    with open(config_filename, 'r') as config_f:
        for line in config_f:
            # skip the comment lines
            if line.strip()[0] == '#':
                continue
            # skip the empty lines
            if line.strip() =='':
                continue
            # parse the parameter lines
            tokens = line.split('=')
            if tokens[0].strip() in par_dict:
                par_dict[tokens[0].strip()] = tokens[1].strip().strip('\'')
            else:
                print("Unknown parameter name: " + tokens[0].strip())
    return par_dict

# Class to change working directory
class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

# Function to print a line of a pdb file
def parse_atom_line(line):
    atom = []
    # ATOM or HETATM
    atom.append(line[0:6].strip())
    # atom serial number
    atom.append(int(line[6:11].strip()))
    # atom name
    atom.append(line[12:16].strip())
    # location indicator
    atom.append(line[16:17].strip())
    # residue name
    #atom.append(line[17:20].strip())
    # All the DUD-e ligands residue names are CHE
    atom.append('CHE')
    # chain ID
    atom.append('')
    # residue sequence number
    #atom.append(line[22:26].strip())
    atom.append(1)
    # code for insertion of residues
    atom.append(line[26:27].strip())
    # coordinate X
    atom.append(float(line[30:38].strip()))
    # coordinate Y
    atom.append(float(line[38:46].strip()))
    # coordinate Z
    atom.append(float(line[46:54].strip()))
    # occupancy
    atom.append(float(line[54:60].strip())) 
    # temperature factor
    atom.append(float(line[60:66].strip())) 
    # element symbol
    atom.append(line[12:16].strip()[0])
    # charge on the atom
    # PDB file doesn't contain charge information
    atom.append('')
    return atom

# Function to write pdb file
def write_pdb(atom, outfile):
    outfile.write("%-6s%5d %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %2s%2s\n" % (atom[0], atom[1],atom[2],atom[3],atom[4],atom[5],atom[6],atom[7],atom[8],atom[9],atom[10],atom[11],atom[12],atom[13],atom[14]))

# Function to modify the 9th line 64th char to 0
def change_prmtop(infilename, outfilename):
    infile = open(infilename, 'r')
    outfile = open(outfilename, 'w')
    for i in range(8):
        outfile.write(infile.readline())
    str1 = ''
    line9list = list(infile.readline())
    line9list[63] = 0
    outfile.write(str1.join(map(str, line9list)))
    while True:
        line = infile.readline()
        if not line: break
        outfile.write(line)
    infile.close()
    outfile.close()


# Main
if __name__ == '__main__':
    main()
