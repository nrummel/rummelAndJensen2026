# Bi-PQN4LCPs

Code to recreate numerical experiments for the [paper](https://arxiv.org/abs/2604.10089) 

# Set Up
In your terminal clone the code and initialize the submodules:
```bash 
git clone git@github.com:nrummel/Bi-PQN4LCPs.git
git submodule update --init --recursive
```

# Recreate Results 
You can recreate our plots by running the following in the specified order. First, in MATLAB, recreate the offline results:
```MATLAB
cd PATH/TO/REPO/src/scripts
allPlots
```
The table of our offline results can be created through Julia:
```julia
cd PATH/TO/REPO/src/scripts
using Pkg; Pkg.activate(".")
juliaPlots
```
Our online results can be plotted through the iPython notebook. This conda environment can be loaded through the following commands in the terminal:
```bash 
cd PATH/TO/REPO/src/scripts`
conda env create -f environment.yml
conda activate Bi-PQN4LCPs
```
Then open the `PATH/TO/REPO/src/scripts/plots.ipynb` and run all the cells with the conda environment specified.


# Running Spherical Collision Code
In order to recreate our results from scratch, you will need to be able to run MATLAB commands on a machine with ~400GB of ram and at least 24 modern CPU's.
## Getting Started
We suggest that others first get a simple test run. Try to run the MATLAB code below first before trying to run larger simulations
```MATLAB
cd PATH/TO/REPO/src/scripts
n                 = 3; % lattice size corresponds to n^3 number of particles
meanRadius        = 1;
Cdst              = 3;
p                 = 8; % Decrease to run on smaller machines, but be warned the accuracy is highly dependent on the discretization of the BIs
ep                = .3;
polydisperseRatio = 0.2;
lcpSlvr           = 'bifi';
lcpTol            = 1e-8;
lcpMaxIter        = 100;
lcpWarmStart      = false;
lcpPLo            = 4;
lcpPLo            = 1e-6; % 1e-8 also works well
Nt                = 200;
dt                = 0.1;
tdisc             = 'euler';
gamma             = 1;
lambda            = 0.1;
mdist             = 3;
denseMV           = true; % switch to false if you want to use FMM
denseforce        = true;
gmresTol          = 1e-8;
saveLCPs          = false;
initMode          = 'lattice';
loadIntermediate  = false;
plotFlag          = false;
seed              = 1;
scaleFlag         = true;
Test_ModLap_Mobility_Amphi(...
    n,meanRadius,Cdst,p,ep,polydisperseRatio, ...
    lcpSlvr,lcpTol,lcpMaxIter,lcpWarmStart, lcpPLo, lcpTolLo, ...
    Nt,dt,tdisc, ...
    gamma,lambda,mdist,denseMV,denseforce,gmresTol, ...
    saveLCPs,initMode,loadIntermediate,plotFlag,seed, ...
    scaleFlag,resDir)
```
## Recrating Our Exact Results
The job scripts to run our code via [slurm](https://slurm.schedmd.com/overview.html) on the [CU Boulder](https://curc.readthedocs.io/en/latest/) Alpine and Blanca clusters are 
- **Offline Results**: Running these files in this order should recreate comparable results to our paper's offline result section.
    - `src/scripts/amphi.lattice.offline.n_5.sh` runs the simulation for the results specified by our offline results section.
    - `src/scripts/saveDenseMatVecs.py` kicks off a family of slurm job arrays to create dense representations of the LCP matrices for each parameterization ($p, \epsilon_\mathrm{gmres}$) at each time step.
    - `src/scripts/zipDenseMats.m` puts all of the dense matrices into one file for convenience.
- **Online Results**: These four scrips will kick off slurm job arrays to recreate results comparable to our papers online results section.
    - `src/scripts/amphi.lattice.online.n_3.sh`
    - `src/scripts/amphi.lattice.online.n_4.sh`
    - `src/scripts/amphi.lattice.online.n_5.sh`
    - `src/scripts/amphi.lattice.online.n_6.sh`