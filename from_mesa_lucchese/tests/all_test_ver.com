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
./mesa_test_${TEST}.job >&mesa_test_${TEST}.o$ver
done

for TEST in 1 2 ; do
./ckohn_test_${TEST}.job >&ckohn_test_${TEST}.o$ver
done


