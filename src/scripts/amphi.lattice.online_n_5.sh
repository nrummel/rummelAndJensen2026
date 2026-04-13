#!/bin/bash
#SBATCH --array=0-2
#SBATCH --nodes=1
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=24
#SBATCH --mem=400G
#SBATCH --account=ucb289_asc4
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --job-name=amphi.lattice.online.n_5
#SBATCH --output=/projects/$USER/onlineResults/amphi.lattice.online.n_5.out
#SBATCH --error=/projects/$USER/onlineResults/amphi.lattice.online.n_5.err 
module purge
## Begin script
echo "=="
echo "||"
echo "|| Begin Execution of fd in slurm batch script."
echo "||" 
echo "=="

## Global Hyper-parameters
meanRadius=1; # mean radius of each of the particles
ep=.3; # distance where we consider collisions in the LCP
p=8; # number of spherical harmonics
Cdst=3.0; # initial distance of the particle centers
polydisperseRatio=0;
# LCP params
lcpSlvr="bbpgd";
lcpTol=0.00000001; # 1e-8
lcpMaxIter=100;
lcpWarmStart=0;
lcpPLo=6;
lcpTolLo=0.000001; # 1e-6
# Time disc params
Nt=100; # number of time-steps
dt=.1; # time discretization
scaleFlag=1; # special initialization to get particles to be close to guarantee the number of collisions to be high near the initial time step

## Define list of parameters
algos=("bbpgd" "proxquasinewton" "bifi")
latticeSize=(5)
polyDisperseRatios=(0.2)
numLatticeSize=1
numAlgos=3

# for SLURM_ARRAY_TASK_ID in {0..17}
# do
## specify ix's from slurm task id
i=$SLURM_ARRAY_TASK_ID
algoIx=$(expr $i % $numAlgos)
ii=$(expr $i / $numAlgos)
latticIx=$(expr $ii % $numLatticeSize)
polyIx=$(expr $ii / $numLatticeSize)
## Get params from ix's
lcpSlvr=$(expr ${algos[$algoIx]})
n=$(expr ${latticeSize[$latticIx]})
polyRatio=$(expr ${polyDisperseRatios[$polyIx]})
echo "lcpSlvr=$lcpSlvr, n=$n, polyRatio=$polyRatio"

module purge
module load matlab gcc 

resDir="/projects/$USER/onlineResults";
# call MATLAB
cd /projects/$USER/Bi-PQN4LCPs/src/scripts
LD_PRELOAD=/curc/sw/install/gcc/14.2.0/lib64/libgfortran.so.5 matlab -nodesktop -nodisplay -r "clear;clc; Test_ModLap_Mobility_Amphi($n,$meanRadius,$Cdst,$p,$ep,$polyRatio,'$lcpSlvr',$lcpTol,$lcpMaxIter,$lcpWarmStart, $lcpPLo, $lcpTolLo,$Nt,$dt,$scaleFlag,'$resDir')); quit;"
echo "- LD_PRELOAD=/curc/sw/install/gcc/14.2.0/lib64/libgfortran.so.5 matlab -nodesktop -nodisplay -r \"clear;clc; Test_ModLap_Mobility_Amphi($n,$meanRadius,$Cdst,$p,$ep,$polyRatio,'$lcpSlvr',$lcpTol,$lcpMaxIter,$lcpWarmStart, $lcpPLo, $lcpTolLo,$Nt,$dt,$scaleFlag,'$resDir'); quit;\""
## Finish scripts
echo "=="
echo "||"
echo "|| Execution of ./main in slurm batch script complete."
echo "||"
echo "=="