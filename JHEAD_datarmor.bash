#!/bin/bash 
XXDEPENDENCYXX
#PBS -l mem=8g
#PBS -l walltime=00:XXTIMEXX:00
#PBS -N XXSPTYPEXX_XXCASEXX_XXSYEAXX_XXEYEAXX
#PBS -e XXINITDIRXX/XXPTYPEXX/JOBS-OUT/job_out_XXCASEXX_XXSYEAXX_XXEYEAXX.e
#PBS -o XXINITDIRXX/XXPTYPEXX/JOBS-OUT/job_out_XXCASEXX_XXSYEAXX_XXEYEAXX.o

cd $PBS_O_WORKDIR

#qstat -f $PBS_JOBID
#echo $TMPDIR
#echo $SCRATCH
#echo $DATAWORK
#echo $HOST
#pbsnodes $HOST

export MODULEPATH="$MODULEPATH:/home1/datahome/ctalandi/modules::.:"

source /usr/share/Modules/3.2.10/init/bash

source /appli/anaconda/versions/miniforge3-24.11.3-0/etc/profile.d/conda.sh
which conda 
conda activate My-JupDask-2025-10

