#!/bin/bash
# put this file run.sh in the directory of your project
# modify MY_PROGRAM_NAME and FILE.m

# bash run.sh
# qstat 
# qnodes
# qdel

LANG=en_US.UTF-8

PNAME=${PWD##*/}
NODE=1
NP=32
QSUBTIME="9600:00:00"
NOWDIR=`pwd`
QSUB=$PNAME.qsub

cat <<EOF >$QSUB
#!/bin/bash

#PBS -q old
#PBS -V 
#PBS -N $PNAME
#PBS -l nodes=$NODE:ppn=$NP
#PBS -l walltime="$QSUBTIME"
#PBS -o $NOWDIR/log.txt
#PBS -j oe 

cd $NOWDIR
/home/software/MATLAB/R2019b/bin/matlab <main.m> record.txt
sleep 1
EOF

qsub $QSUB

