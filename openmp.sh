#! /bin/sh
g++ $1 -fopenmp 
filename="$1.sh"
mv job.sh filename
qsub -g tga-hpc-lecture filename
mv filename job.sh
