#! /bin/sh
g++ $1 -fopenmp -o openmp
filename="$2.sh"
mv job.sh $filename
qsub -g tga-hpc-lecture $filename
mv $filename job.sh
