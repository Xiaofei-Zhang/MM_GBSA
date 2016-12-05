# get_avg_energy_best_5.py
# get the Energy of the MMPBSA results
# Usage
#	python get_avg_energy_best_5.py $LIG_LIST $RCPTR
#
# $LIG_LIST is the file contains all the ligand names
# $RCPTR is the receptor name
#
# Author: Xiaofei Zhang
# Date: July 25 2016

import os
import sys

def main():
        if len(sys.argv) != 3:
		print("Usage: python get_energy.py $LIG_LIST $RCPTR")
		sys.exit(-1)
	lig_list = [line.strip() for line in open(sys.argv[1], 'r')]
	rcptr = sys.argv[2]
	# DIR contains all the folders
	job_dir = '/scratch/xzh289/MM_PBSA_data/1QCF_20/out'
	model_number = 20
	avg_number = 5

	energy_of_pdbqt = {} 
	file_names = []
	for lig in lig_list:
		energy_of_models = []
		base_name = rcptr + '_' + lig
		pdbqt_name = base_name + '.pdbqt'
		file_names.append(pdbqt_name)
		for i in range(1, model_number + 1):
			dir_name = job_dir + '/' + base_name + '.' + str(i)
			if not os.path.isfile(dir_name + '/FINAL_RESULTS_MMPBSA.dat'):
				continue
			energy_of_models.append(get_one_energy(dir_name))
		if len(energy_of_models) == 0:
                        energy_of_pdbqt[pdbqt_name]= 999999.99
			continue
		energy_of_pdbqt[pdbqt_name] = avg_list(energy_of_models, avg_number)
	
	out_fname = rcptr + '_mmpbsa_energy_avg_best_5'
	with open(out_fname, 'w') as out_file:
		out_file.write('Order\tEnergy\tLigand')
		order = 1
		for pdbqt in file_names:
			out_file.write('\n'+str(order)+'\t' +str(energy_of_pdbqt[pdbqt])+ '\t' + pdbqt)
			order = order + 1
	
	



def avg_list(mylist, models):
	sum_val = 0
	mylist.sort(reverse=True)
	nums = 0
	if models > len(mylist):
		nums = len(mylist)
	else:
		nums = models
	for i in range(nums):
		sum_val = sum_val + mylist[i]
	avg_val = sum_val / nums
	return avg_val

def get_one_energy(dir_name):
	infile = open(dir_name + '/FINAL_RESULTS_MMPBSA.dat', 'r')
	energy = 0
	for line in infile:
		if line[0:11] == 'DELTA TOTAL':
			energy = float(line.split()[2])
			break
	infile.close()
	return energy


# Main
if __name__ == '__main__':
	main()
