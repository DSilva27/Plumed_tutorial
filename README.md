# Plumed_tutorial

Scrips to set-up a multi-replica plumed MD simulation using the RMSD value as a collective variable. Installation parameters and codes are based on [plumed_em_md from fraser-lab](https://github.com/fraser-lab/plumed_em_md).

## Requires:

* python 3.5 or higher
* PLUMED v2.6
* Gromacs 2018.8

## Install

### Installing PLUMED

Note: do this somewhere else!
```bash
  git clone https://github.com/plumed/plumed2.git
  cd plumed2
  git checkout v2.6
  ./configure --disable-python
  make -j 12
```

### Installing GROMACS
Note: do this somewhere else!
```bash
  wget http://ftp.gromacs.org/pub/gromacs/gromacs-2018.8.tar.gz
  tar -xvf gromacs-2018.8.tar.gz
  cd gromacs-2018.8
  plumed-patch -p --shared
  ## Choose gromacs-2018.8
  
  mkdir build
  cd build
  mkdir [Install folder]/gromacs-2018.8-bin/
  
  # These parameters are still to be defined, I'll leave the ones from fraser-lab for now.
  cmake ../ -DBUILD_SHARED_LIBS=ON -DGMX_OPENMP=OFF -DGMX_THREAD_MPI=OFF -DGMX_GPU=OFF 
  -DCMAKE_INSTALL_PREFIX=[Install folder]/gromacs-2018.8-bin/ -DCMAKE_CXX_COMPILER=mpic++ 
  -DCMAKE_C_COMPILER=mpicc -DGMX_MPI=ON -DGMX_USE_RDTSCP=off
  
  make -j 12
  make install
```

### Configuring .bashrc
```bash
  # Add to .bashrc
  source ~/plumed2/sourceme.sh
  source ~/gromacs-2016.5-bin/bin/GMXRC
```

## Prepare the simulation:

## Set-up simulation folder
```bash
  mkdir simulation
  cp MDP/* simulation
  cp python_files/* simulation
  cp structural_files/* simulation
  cd simulation
  ```

### Prepare production MD

Notes: 
* Aside from min_and_equib.py, the other codes are really quick. For prep_plumed.py one core should suffice.
* If you are using the pdb files used in this repo, then `input_pdb` should be `apo.pdb` and `ref_pdb` should be `holo.pdb`.

```bash
  python min_and_equib.py --input_pdb [path to input pdb] --n_proc [number of processors] --n_rep [number of replicas]
  python prep_replicas.py --n_rep [number of replicas]
  python prep_plumed.py --ref_pdb [pdb used for RMSD calculation] --n_proc [number of processors] --n_rep [number of replicas]
  ```
  
### Run production MD

```bash
  mpirun -n [n_ranks] gmx_mpi mdrun -plumed -multi [n_rep] -ntomp [n_threads]
  ```
  
### Restarting
Uncomment out #RESTART in plumed.dat
```bash
  mpirun -n [n_ranks] gmx_mpi mdrun -plumed -multi [n_rep] -ntomp [n_threads] -cpi state
  ```
