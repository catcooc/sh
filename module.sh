#BSUB -q e52692v2ib
#BSUB -n 48
#BSUB -R largemem
#BSUB -J MgSiO3
#BSUB -o out
#BSUB -e err
module load ips/2018u4
mpiexec ./dqmc.x
