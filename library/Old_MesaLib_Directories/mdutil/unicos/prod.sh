#!/bin/csh
onintr wrapup
set cfsdir = /112793
set jobid = $QSUB_REQNAME
ja
date > $jobid.acct
mkdir t
setenv MESA_TMP t
setenv MESA_BIN /usr/tmp/russo/bin
cfs get $MESA_BIN/mesa:$cfsdir/mesaexec/mesa
cfs get Cu.d10s

mesa inp=Cu.d10s out=Cu.d10s.out chk=mesachk saveint=t/foo >>& mesa.log

cfs store mesa.log
cfs store Cu.d10s.out
cfs store mesachk:Cu.d10s.chk
cfs store t/foo:Cu.ints

wrapup:
ja -cfst >> $jobid.acct

exit

