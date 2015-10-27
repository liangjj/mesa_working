#!/bin/bash

prg=
ver=
cmnd=
if [[ $1 != "" ]] ; then
case $1 in
     all | lib | clean )
           cmnd=$1
           ver=$2
           ;;
     .* )
           cmnd=all
           ver=$1
           ;;
     one )
           cmnd=one 
           prg=$2
           ver=$3
           ;;
     tar )
           cmnd=tar 
           ;;
     tarc )
           cmnd=tarc 
           ;;
     * )
           cmnd=one
           prg=$1
           ver=$2
           ;;
esac
else
   cmnd=all
   ver=
   prg=
fi
echo command is $cmnd 
echo program is $prg 
echo version is $ver

if [[ $cmnd == tar ]] ; then
  if [[ $MESA == "" ]] ; then
     MESA=$PWD
  fi
  MESAMAIN=${MESA##*/}
  cd $MESA
  cd ..
  tar -c -z -f ${MESAMAIN}.`date +%d.%m.%Y`.tgz --exclude "$MESAMAIN/bin*" \
  --exclude "*.a" --exclude "*/map.out" --exclude "*.mod" --exclude "*.o" \
  $MESAMAIN 
elif [[ $cmnd == tarc ]] ; then
  if [[ $MESA == "" ]] ; then
     MESA=$PWD
  fi
  MESAMAIN=${MESA##*/}
  cd $MESA
  cd ..
  tar -c -z -f ${MESAMAIN}.`date +%d.%m.%Y`c.tgz --exclude "$MESAMAIN/bin*" \
  --exclude "*.a" --exclude "*.v*" --exclude "*/map.out" --exclude "*~" --exclude "*.mod" \
  --exclude "*.o" --exclude "*.o.*" $MESAMAIN 
elif [[ $cmnd == clean ]] ; then
   cd library
   ./Makemesalib.sh clean $ver
   cd ..
   cd dir_mylib
   rm -f mylib.a
   if [[ $ver != "" ]] ; then
      rm -f mylib$ver.a
   else
      rm -f mylib*.a
   fi
   cd ..
   rm -r -f bin
   if [[ $ver != "" ]] ; then
      rm -r -f bin$ver
   else
      rm -r -f bin*
   fi
else
if [[ $cmnd == library || $cmnd == all ]] ; then
rm -f -r bin$ver
mkdir bin$ver
fi

if [[ $ver != "" ]] ; then
  rm -r -f bin
  ln -s bin$ver bin
fi

source include/$MACH$ver.sh

cp optmesa bin

cd library

if [[ $cmnd == library || $cmnd == all ]] ; then
echo ' Making library'
rm -f *.o
./Makemesalib.sh make $ver
rm -f *.o
fi
./Makemesalib.sh connect $ver
cd ..

#
cd dir_mylib
if [[ $cmnd == library || $cmnd == all ]] ; then
   echo ' Making mylib'
   rm -f *.o
   rm -f *.a
   make
   rm -f *.o
   if [[ $ver != "" ]] ; then
      mv mylib.a mylib$ver.a
   fi
fi

if [[ $ver != "" ]] ; then
   rm -f mylib.a
   ln -s mylib$ver.a mylib.a
fi
cd ..

if [[ $cmnd == all ]] ; then
   for PROG in m[0-9]* k* ; do
     cd $PROG
     echo " Making $PROG"
     rm -f *.o
     make
     rm -f *.o
     cd ..
   done
else
  PROG=$prg
  cd $PROG
  echo " Making $PROG"
  rm -f *.o
  make all
  rm -f *.o
  cd ..
fi

if [[ $ver != "" ]] ; then
  rm -r -f bin
  rm dir_mylib/mylib.a
  cd library
  ./Makemesalib.sh disconnect
  cd ..
fi
fi












