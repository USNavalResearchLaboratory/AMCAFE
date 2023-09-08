#! /bin/bash
# Request maximum wallclock time for the job (HH:MM:SS)

# Job size: select=# of nodes, ncpus=cores/node (this is system), mpiprocs=MPI_procs/node (this is what you want)
# #PBS -l select=4:ncpus=8:mpiprocs=8+1:ncpus=1:mpiprocs=1

# # THIS IS PRIME TIME RUN
#PBS -l walltime=10:00:00
#PBS -l select=2:ncpus=48:mpiprocs=48
#PBS -q standard

# # THIS IS DEBUG RUN
# #PBS -l walltime=00:30:00
# #PBS -l select=1:ncpus=48:mpiprocs=48
# #PBS -q debug

# Specify job name
#PBS -N sdx1a
# Select Project ID (from pIE allocation)                                                                                                                         
#PBS -A NRLDC04161571
# Pass environment variables to the job
#PBS -V
# Join output files
#PBS -j oe                                                                                                                                                      
module swap compiler/intel/2019.4.243 compiler/intel/2019.5.281
module swap mpt/2.20 compiler/intelmpi/2019.5.281
module load hdf5-parallel/intel-18.1.163/1.10.5
module load gcc/9.2.0
cd /p/work1/kteferra
cp /p/home/kteferra/Documents/CA/AMCAFE2/AMCAFE/scripts/SDX1a_.in /p/work1/kteferra
# must make sure # of processors is consistent with what's being requested above
/p/app/intel/parallel_studio_xe_2019_update5/impi/2019.5.281/intel64/bin/mpirun -n 96 /p/home/kteferra/Documents/CA/AMCAFE2/AMCAFE/src/src-mpi/cafe SDX1a_.in


