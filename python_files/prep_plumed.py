import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--map_dat", help="full path to the GMM dat file for the map",required=True)
parser.add_argument("--n_proc", help="number of threads to use",required=True)
parser.add_argument("--n_rep", help="number of replicas to simulate",required=True)
parser.add_argument("--res", help="Resolution of the cryo-em map", required=True, type=float)
parser.add_argument("--sigma_min", help="SIGMA_MIN parameter for EMMI", default=0.02)

args = parser.parse_args()                                                                                                                                          
map_plumed_dat_file = args.map_dat                                                                                                                                  
n_proc = int(args.n_proc)                                                                                                                                           
replicas = int(args.n_rep)  
res = args.res
sigma_min = args.sigma_min

import os

try:
    from subprocess import STDOUT, check_output, CalledProcessError
except ImportError:  # pragma: no cover
    # python 2.6 doesn't include check_output
    # monkey patch it in!
    # from: https://stackoverflow.com/questions/4814970/subprocess-check-output-doesnt-seem-to-exist-python-2-6-5
    import subprocess
    STDOUT = subprocess.STDOUT

    def check_output(*popenargs, **kwargs):
        if 'stdout' in kwargs:  # pragma: no cover
            raise ValueError('stdout argument not allowed, '
                             'it will be overridden.')
        process = subprocess.Popen(stdout=subprocess.PIPE,
                                   *popenargs, **kwargs)
        output, _ = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            raise subprocess.CalledProcessError(retcode, cmd,
                                                output=output)
        return output
    subprocess.check_output = check_output

    # overwrite CalledProcessError due to `output`
    # keyword not being available (in 2.6)
    class CalledProcessError(Exception):

        def __init__(self, returncode, cmd, output=None):
            self.returncode = returncode
            self.cmd = cmd
            self.output = output

        def __str__(self):
            return "Command '%s' returned non-zero exit status %d" % (
                self.cmd, self.returncode)
    subprocess.CalledProcessError = CalledProcessError

out = str(check_output(["sed '3q;d' conf_box_oriented.gro"],shell=True))
ref = ",".join(out.split()[-3:])

#create protein-no-negative group in index.ndx
conf_in = open("conf_box_oriented.gro")
negative_charges = ["OD1","OD2","OE1","OE2"]
negative_residues = ["ASP","GLU"]
negative_atoms = []
for line in conf_in:
    try:
        if line.split()[1] in negative_charges and line.split()[0][-3:] in negative_residues:
            negative_atoms.append(line.split()[2])
    except:
        pass

os.system("cp index.ndx index.ndx_BACKUP")
index_file = open("index.ndx","r+")
in_protein = 0
protein_atoms = []
for line in index_file:
    if "[ C-alpha ]" in line:
        in_protein = 0
    if in_protein:
        protein_atoms.extend(line.split())
    if "[ Protein-H ]" in line:
        in_protein = 1

for negative_atom in negative_atoms:
    protein_atoms.remove(negative_atom)

index_file.write("[ Protein-no-negative-no-h ]\n")
i = 1
newline = ""
for atom in protein_atoms:
    newline = newline + atom + " "
    i +=1
    if i % 15 == 0:
        index_file.write(newline)
        index_file.write("\n")
        newline = ""
index_file.write(newline)
index_file.close()

from multiprocessing import Pool

def create_topols(i):
    os.system("gmx_mpi grompp -f nvt_2016.mdp -c conf_{i}.gro -o topol{i}.tpr -maxwarn 3".format(i=i))

p = Pool(n_proc)
print(p.map(create_topols,range(0,replicas,1)))

fout = open("plumed.dat","w")
fout.write("""# RESTART

# include topology info: this is needed to identify atom types
MOLINFO STRUCTURE=structure.pdb

# define all heavy atoms using GROMACS index file
# which can be created with gmx_mpi make_ndx
protein-h: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H
protein: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein
protein-no-negative: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-no-negative-no-h
# water: GROUP NDX_FILE=index.ndx NDX_GROUP=Water

# make protein whole: add reference position of first heavy atom (in nm)
WHOLEMOLECULES ADDREFERENCE ENTITY0=protein REF0={ref}

# create EMMI score
EMMI ...
LABEL=gmm NOPBC TEMP=300.0 NL_STRIDE=50 NL_CUTOFF=0.01
ATOMS=protein-no-negative GMM_FILE={map_dat}
SIGMA_MIN={sigma_min} RESOLUTION={res} NOISETYPE=MARGINAL
...

3# translate into bias                                                                                                                                                
34emr: BIASVALUE ARG=gmm.scoreb STRIDE=1                                                                                                                               
35                                                                                                                                                                     
36SELECTOR NAME=BETA VALUE=0                                                                                                                                           
37                                                                                                                                                                     
38RESCALE ...                                                                                                                                                          
39LABEL=res ARG=emr.bias TEMP=300                                                                                                                                      
40SELECTOR=BETA MAX_RESCALE=1000 NBIN=20                                                                                                                               
41W0=1000 BIASFACTOR=150000 BSTRIDE=1000 BFILE=bias.dat                                                                                                                
42MC_STRIDE=1                                                                                                                                                          
43...     

# print useful info to file
PRINT ARG=* FILE=COLVAR STRIDE=1000
""".format(ref=ref, map_dat=map_plumed_dat_file, sigma_min=sigma_min, res=res))
fout.close()

print("""
You should be set to run:
mpirun -n N_RANKS gmx_mpi mdrun -plumed -multi {n_rep} -ntomp N_THREADS
if you need to restart:
Edit plumed.dat
Uncomment out #RESTART
Then on command line try:
mpirun -n N_RANKS gmx_mpi mdrun -plumed -multi {n_rep} -ntomp N_THREADS -cpi state
you can analyze on the fly or at end with:
generate_trajectories_as_pdbs.py
""".format(n_rep=replicas))

