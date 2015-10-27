
# %W% %G%
#  This shell script runs several Mesa input decks one after the other.
#  To use this shell, first prepare an input file with the names of
#  the input decks to be run.  It might look like:
#
#  test1
#  test2
#  end
#
#  Suppose you call it test.inp; the shell is then invoked by the command
#  prodmesa < tests.inp
#
#  the following variables should be set to your own personal run environment
rundir=/usr2/rlm/run
inpdir=/usr2/rlm/inp
outdir=/usr2/rlm/out
chkdir=/usr2/rlm/chk
#
#
   cd $rundir
   read name
   until test $name = end
      do
         echo executing $name
         mesa inp=$inpdir/$name out=$outdir/$name 
         read name
      done

