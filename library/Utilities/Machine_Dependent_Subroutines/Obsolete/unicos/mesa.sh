
# 	%W%   %G%
#
#    	script for running the mesa links.
#       the first line of the script is deliberately left blank in order
#       to invoke the Bourne shell.
#
#       set up a temporary directory n the directory where this will run 
#       called "t".  the path is set up to point to wherever the executbles
#       reside, and the tmpdir variable set to a large temporary space
#       for scratch files.  note that you may wish to modify the subroutine
#       lxopen.f in order to set up the default file names as you wish.
#       otherwise, they can all be defaulted here and passed in the
#       execution line.
#
#
PATH=/usr/tmp/mob/mesa/bin:/usr/bin:/usr/ucb:/bin:/usr/local/bin
export PATH
TMPDIR=/usr/tmp
export TMPDIR
#
#       set up the default file names.
inp='mesa.inp'
out='mesa.out'
chk='t/chk'
dat='/usr/tmp/mesa/mesa.dat'
int='t/int'
rwf='t/rwf'
siz='1500000'
rint='t/rint'
tint='t/tint'
gint='t/gint'
rdint='t/rdint'
dint='t/dint'
zint='t/zint'
ham='t/ham'
moden='t/moden'
aoden='t/aoden'
saoden='t/saoden'
gden='t/gden'
fci='t/fci'
#
kohn='t/kohn'
kohndt='t/kohndt'
grid='t/grid'
orbs='t/orbs'
vstat='t/vstat'
ylms='t/ylms'
bessel='t/bessel'
knints='t/knints'
tmat='t/tmat'
blktmt='t/blktmt'
#
#       set up the default initial and terminal links.
start='m0'
stop='m998'
#
#       check to see if any of the files have been replaced on the
#       command line.
for file in $*
do
   case $file in
      inp=*) 
        inp=`expr $file : 'inp=\(.*\)'` ;;
      out=*) 
        out=`expr $file : 'out=\(.*\)'` ;;
      chk=*) 
        chk=`expr $file : 'chk=\(.*\)'` ;;
      dat=*) 
        dat=`expr $file : 'dat=\(.*\)'` ;;
      rwf=*) 
        rwf=`expr $file : 'rwf=\(.*\)'` ;;
      int=*) 
        int=`expr $file : 'int=\(.*\)'` ;;
      rint=*) 
        rint=`expr $file : 'rint=\(.*\)'` ;;
      tint=*) 
        tint=`expr $file : 'tint=\(.*\)'` ;;
      gint=*) 
        gint=`expr $file : 'gint=\(.*\)'` ;;
      rdint=*) 
        rdint=`expr $file : 'rdint=\(.*\)'` ;;
      dint=*) 
        dint=`expr $file : 'dint=\(.*\)'` ;;
      zint=*) 
        zint=`expr $file : 'zint=\(.*\)'` ;;
      ham=*) 
        ham=`expr $file : 'ham=\(.*\)'` ;;
      moden=*) 
        moden=`expr $file : 'moden=\(.*\)'` ;;
      aoden=*) 
        aoden=`expr $file : 'aoden=\(.*\)'` ;;
      saoden=*) 
        saoden=`expr $file : 'saoden=\(.*\)'` ;;
      gden=*) 
        gden=`expr $file : 'gden=\(.*\)'` ;;
      fci=*) 
        fci=`expr $file : 'fci=\(.*\)'` ;;
      start=*) 
        start=`expr $file : 'start=\(.*\)'` ;;
      stop=*) 
        stop=`expr $file : 'stop=\(.*\)'` ;;
   esac
done
#
#       run the links. 
next=$start
#
#       run until the next message is not a link number 
until echo $next | egrep -v 'm[0-9]+'
do
#       is this the terminal link?
        if echo $next | egrep $stop
        then
                break
        fi
#
        echo $next inp=$inp out=$out chk=$chk rwf=$rwf int=$int dat=$dat siz=$siz 
        next=`$next inp=$inp out=$out chk=$chk rwf=$rwf int=$int dat=$dat  siz=$siz`
done
#
#       we have finished.  
