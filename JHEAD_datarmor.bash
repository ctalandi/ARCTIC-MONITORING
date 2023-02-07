#!/bin/bash 
XXDEPENDENCYXX
#PBS -l mem=60g
#PBS -l walltime=00:XXTIMEXX:00
#PBS -N XXSPTYPEXX_XXCASEXX_XXSYEAXX_XXEYEAXX
#PBS -e XXINITDIRXX/XXPTYPEXX/JOBS-OUT/job_out_XXSYEAXX_XXEYEAXX.e
#PBS -o XXINITDIRXX/XXPTYPEXX/JOBS-OUT/job_out_XXSYEAXX_XXEYEAXX.o

cd $PBS_O_WORKDIR

#qstat -f $PBS_JOBID
#echo $TMPDIR
#echo $SCRATCH
#echo $DATAWORK
#echo $HOST
#pbsnodes $HOST

source /usr/share/Modules/3.2.10/init/bash
module purge
module load   NETCDF/4.3.3.1-mpt-intel2016
module load   nco/4.6.4_gcc-6.3.0
module load   vacumm/3.4.0
module list

