#!/bin/bash
#PBS -l nodes=2:ppn=4

echo "------------------"
echo

if [ -n "$PBS_O_WORKDIR" ]
then
   cd $PBS_O_WORKDIR
   echo "Job working directory: $PBS_O_WORKDIR"
else
   echo "Job working directory: "`pwd`
fi
echo

if [ -n "$PBS_NODEFILE" ]
then
   num=`cat $PBS_NODEFILE | wc -l`
   echo "Requested processors: $num"
   echo "Nodes:"
   uniq $PBS_NODEFILE
   echo
else
   echo "No node file - running locally"
fi

echo "Job starting at `date`"
echo

START=`date '+%s'`

# The $val variable is passed to this script with qsub's -v option, e.g.,
#   % qsub -v val=2 run_code4.m

matlab -nodisplay  -nodesktop -nosplash -r "evscript"
# matlab -nodisplay -nojvm -nodesktop -nosplash -singleCompThread -r "makemode"

END=`date '+%s'`

echo
echo "Job finished at `date`"
echo
echo 'Total Execution time: '`expr $END - $START`' seconds'
echo

