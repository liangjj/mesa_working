#!/bin/bash

if [ $MESA == "" ] ; then
   echo Please define MESA
   exit
fi
cd $MESA
rm other/MoreThan72.out
for FILE in $(find . -name "*.f" -print | grep -v "\.v" | grep -v "./other" ) ; do
echo File $FILE >> other/MoreThan72.out
other/MoreThan72.exe <$FILE >>other/MoreThan72.out
done
for FILE in $(find . -name "*.for" -print | grep -v "\.v" | grep -v "./other" ) ; do
echo File $FILE >> other/MoreThan72.out
other/MoreThan72.exe <$FILE >>other/MoreThan72.out
done

