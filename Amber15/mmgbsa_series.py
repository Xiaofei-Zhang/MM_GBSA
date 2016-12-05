import os,sys


def main():

    out_dir = '/oasis/scratch/comet/xzh33/temp_project/MM_GBSA_Benchmarking/1N/out/'
    cmd = 'MMPBSA.py -O -i ~/MMGBSA_benchmarking/mmpbsa_in/mmpbsa.in -o FINAL_RESULTS_MMPBSA.dat -cp cmplx_i.prmtop -rp rcptr_i.prmtop -lp lig_i.prmtop -y cmplx.dcd'
    lig_list = [line.strip() for line in open(sys.argv[1], 'r')]
    rcptr_list = [line.strip() for line in open(sys.argv[2], 'r')]
    job_dict = []

    for rcptr in rcptr_list:
        for lig in lig_list:
            job_dict.append(rcptr + '_' + lig)

    dirs = os.listdir(out_dir)
    for mydir in dirs:
        if mydir.split('.')[0] in job_dict:
            with cd(out_dir+'/'+mydir):
                os.system(cmd)


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


# Main
if __name__ == '__main__':
    main()
