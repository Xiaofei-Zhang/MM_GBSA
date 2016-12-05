# change_rcptri.py
# Script to change the receptor pdb file (Used in docking) to amber format pdb 
# that can be used in the MM/GBSA analysis
# Usage:
#         python change_rcptri.py $PDB
#
# $PDB is the receptor pdb file that used in docking
# This script is specificly used to modify the receptor pdb files of the 4 LCK 
# models(1QCF,1Y57,2PL0,2SRC)
# Author: Xiaofei Zhang
# Date: Nov 25 2016



import os,sys

def main():
    if len(sys.argv) != 2:
        print('Usage: python change_rcptri.py $PDB')
        sys.exit(-1)
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        print('File not exist: '+ infilename)
        sys.exit(-1)
    
    basename = infilename[0:-4]
    tmpfilename = basename+ '_tmp.pdb'
    tmpfilename2 = basename+ '_tmp2.pdb'
    outfilename = basename + '_out.pdb'
    # load the amber module on DLX
    cmd1 = 'module load amber/openmpi/12p12'
    os.system(cmd1)
    print cmd1
    # trim the hydrogen atoms in the input pdb file
    cmd2 = 'reduce -Trim ' + infilename + ' > ' + tmpfilename
    os.system(cmd2)
    
    # find out the residue id of phosphotyrosine (Y2P)
    str_Y2P_idx = 'XXXX'
    with open(infilename, 'r') as infile:
        for line in infile:
            if line[13:20] == 'P1  TYR':
                str_Y2P_idx = line[20:26]
                break
    print 'The residue id if phosphotyrosine is :'
    print str_Y2P_idx

    y2p_str=[]
    with open(tmpfilename, 'r') as infile:
        with open(outfilename, 'w') as outfile:
            for line in infile:
                # change the atom name CD to CDI in ILE residue
                if line[13:20] == 'CD  ILE':
                    newline = line[0:13] + 'CD1 ILE' + line[20:]
                    outfile.write(newline)
                    continue
                # change the residuname HSE to HID
                if line[17:20] == 'HSE':
                    newline = line[0:17] + 'HID' + line[20:]
                    outfile.write(newline)
                    continue
                # change the atom name OT1 to O in PRO residue
                if line[13:20] == 'OT1 PRO':
                    newline = line[0:13] + 'O   PRO' + line[20:]
                    outfile.write(newline)
                    continue
                # change the atom name OT2 to OXT in PRO residue
                if line[13:20] == 'OT2 PRO':
                    newline = line[0:13] + 'OXT PRO' + line[20:]
                    outfile.write(newline)
                    continue
                # change the residue name of phosphotyrosine to Y2P, modify the atom name accordingly
                if line[17:26] == 'TYR' + str_Y2P_idx:
                    newline = line[0:17] + 'Y2P' + line[20:]
                    atom_name = newline[13:16]
                    if atom_name == 'OH ':
                        newline = newline[0:13] + 'OG ' + newline[16:]
                    elif atom_name == 'P1 ':
                        newline = newline[0:13] + 'P  ' + newline[16:]
                    elif atom_name == 'O2 ':
                        newline = newline[0:13] + 'O1P' + newline[16:]
                    elif atom_name == 'O3 ':
                        newline = newline[0:13] + 'O2P' + newline[16:]
                    elif atom_name == 'O4 ':
                        newline = newline[0:13] + 'O3P' + newline[16:]
                    y2p_str.append(newline)
                    if atom_name == 'O  ':
                        for i in range(6):
                            outfile.write(y2p_str[i])
                        outfile.write(y2p_str[12])
                        outfile.write(y2p_str[13])
                        outfile.write(y2p_str[6])
                        outfile.write(y2p_str[7])
                        outfile.write(y2p_str[8])
                        outfile.write(y2p_str[9])
                        outfile.write(y2p_str[10])
                        outfile.write(y2p_str[11])
                        outfile.write(y2p_str[14])
                        outfile.write(y2p_str[15])
                    continue
                outfile.write(line)


    '''
    # change CD to CD1 in ILE residues
    cmd3 = 'sed -i \'s/CD  ILE/CD1 ILE/g\' ' + outfilename
    os.system(cmd3)
    # change HSE to HID
    cmd4 = 'sed -i \'s/HSE U/HID U/g\' ' + outfilename
    os.system(cmd4)
    # change OT1 in PRO residue to O
    cmd5 = 'sed -i \'s/OT1 PRO/O   PRO/g\' ' + outfilename
    os.system(cmd5)
    # change OT2 in PRO residue to OXT
    cmd6 = 'sed -i \'s/OT2 PRO/OXT PRO/g\' ' + outfilename
    os.system(cmd6)

    str_Y2P_idx = 'XXXX'
    with open(infilename, 'r') as infile:
        for line in infile:
            if line[13:20] == 'P1  TYR':
                str_Y2P_idx = line[20:26]
                break
    print str_Y2P_idx
    
    cmd7 = 'sed -i \'s/TYR'+ str_Y2P_idx +'/Y2P' + str_Y2P_idx+'/g\' ' + outfilename
    os.system(cmd7)
'''

#    os.system('rm ' + tmpfilename)
if __name__ == '__main__':
    main()


