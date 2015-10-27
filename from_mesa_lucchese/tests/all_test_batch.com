#!/bin/bash
#
if [ "$MESA" == "" ]; then
   export MESA=${PWD/\/tests/}
fi

export ver=$1
if [ "$ver" == "" ] ; then
   echo define version
   exit
fi
echo run tests using version $ver
for TEST in 1 2 3 4 5 6 7; do
qsub -v test=$TEST,ver=$ver,MESA=$MESA mesa_test_batch$MACH.job 
done

for TEST in 1 2 ; do
qsub -v test=$TEST,ver=$ver,MESA=$MESA ckohn_test_batch$MACH.job 
done
