#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=100
#SBATCH --mem=100G
#SBATCH --time=7-00:00:00
#SBATCH --account=blanca-becker
#SBATCH --qos=blanca-becker                          
#SBATCH --job-name=amphi.lattice.offline.n_5
#SBATCH --output=/projects/$USER/offlineResults/amphi.lattice.offline.n_5.out
#SBATCH --error=/projects/$USER/offlineResults/amphi.lattice.offline.n_5.err 

module purge
module load matlab
# These parameters are the main ones that change the scenario
n=5; # size of the lattice
p=8; # number of spherical harmonics
Cdst=2.5; # initial distance of the particle centers
# These parameters should stay the same for the most part
lambda=0.1;
rd=1; # radius of each of the particles
ep=.3; # distance where we consider collisions in the LCP
Nt=250; # number of time-steps
dt=.1; # time discretization
tdisc="euler"; # the only working option I believe, there is some code for AB method but I think it is un tested
saveLCPs=1; # Save the components of each LCP solve
meanRadius=1; # mean radius of each of the particles
polyRatio=0;

# LCP params
lcpSlvr="proxquasinewton";
lcpTol=0.00000001; # 1e-8
lcpMaxIter=100;
lcpWarmStart=0;
lcpPLo=6;
lcpTolLo=0.000001; # 1e-6
# Time disc params
Nt=100; # number of time-steps
dt=.1; # time discretization
resDir="/projects/$USER/offlineResults";

# call MATLAB
cd /projects/$USER/Bi-PQN4LCPs/src/scripts
LD_PRELOAD=/curc/sw/install/gcc/14.2.0/lib64/libgfortran.so.5 matlab -nodesktop -nodisplay -r "clear;clc; Test_ModLap_Mobility_Amphi($n,$meanRadius,$Cdst,$p,$ep,$polyRatio,'$lcpSlvr',$lcpTol,$lcpMaxIter,$lcpWarmStart, $lcpPLo, $lcpTolLo,$Nt,$dt,'$resDir')); quit;"
echo "- LD_PRELOAD=/curc/sw/install/gcc/14.2.0/lib64/libgfortran.so.5 matlab -nodesktop -nodisplay -r \"clear;clc; Test_ModLap_Mobility_Amphi($n,$meanRadius,$Cdst,$p,$ep,$polyRatio,'$lcpSlvr',$lcpTol,$lcpMaxIter,$lcpWarmStart, $lcpPLo, $lcpTolLo,$Nt,$dt,'$resDir'); quit;\""
# done

## Finish scripts
echo "=="
echo "||"
echo "|| Execution of ./main in slurm batch script complete."
echo "||"
echo "=="