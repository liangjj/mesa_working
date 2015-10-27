*deck dchkw
      subroutine dchkw (name, lociw, leniw, locw, lenw, ierr, iter, err)
c***begin prologue  dchkw
c***subsidiary
c***purpose  slap work/iwork array bounds checker.
c            this routine checks the work array lengths and interfaces
c            to the slatec error handler if a problem is found.
c***library   slatec (slap)
c***category  r2
c***type      double precision (schkw-s, dchkw-d)
c***keywords  error checking, slap, workspace checking
c***author  seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c             seager@llnl.gov
c***description
c
c *usage:
c     character*(*) name
c     integer lociw, leniw, locw, lenw, ierr, iter
c     double precision err
c
c     call dchkw( name, lociw, leniw, locw, lenw, ierr, iter, err )
c
c *arguments:
c name   :in       character*(*).
c         name of the calling routine.  this is used in the output
c         message, if an error is detected.
c lociw  :in       integer.
c         location of the first free element in the integer workspace
c         array.
c leniw  :in       integer.
c         length of the integer workspace array.
c locw   :in       integer.
c         location of the first free element in the double precision
c         workspace array.
c lenrw  :in       integer.
c         length of the double precision workspace array.
c ierr   :out      integer.
c         return error flag.
c               ierr = 0 => all went well.
c               ierr = 1 => insufficient storage allocated for
c                           work or iwork.
c iter   :out      integer.
c         set to zero on return.
c err    :out      double precision.
c         set to the smallest positive magnitude if all went well.
c         set to a very large number if an error is detected.
c
c***references  (none)
c***routines called  d1mach, xermsg
c***revision history  (yymmdd)
c   880225  date written
c   881213  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   900805  changed xerrwv calls to calls to xermsg.  (rwc)
c   910411  prologue converted to version 4.0 format.  (bab)
c   910502  corrected xermsg calls to satisfy section 6.2.2 of ansi
c           x3.9-1978.  (fnf)
c   910506  made subsidiary.  (fnf)
c   920511  added complete declaration section.  (wrb)
c   921015  added code to initialize iter and err when ierr=0.  (fnf)
c***end prologue  dchkw
c     .. scalar arguments ..
      double precision err
      integer ierr, iter, leniw, lenw, lociw, locw
      character name*(*)
c     .. local scalars ..
      character xern1*8, xern2*8, xernam*8
c     .. external functions ..
      double precision d1mach
      external d1mach
c     .. external subroutines ..
      external xermsg
c***first executable statement  dchkw
c
c         check the integer workspace situation.
c
      ierr = 0
      iter = 0
      err = d1mach(1)
      if( lociw.gt.leniw ) then
         ierr = 1
         err = d1mach(2)
         xernam = name
         write (xern1, '(i8)') lociw
         write (xern2, '(i8)') leniw
         call xermsg ('slatec', 'dchkw',
     $      'in ' // xernam // ', integer work array too short.  ' //
     $      'iwork needs ' // xern1 // '; have allocated ' // xern2,
     $      1, 1)
      endif
c
c         check the double precision workspace situation.
      if( locw.gt.lenw ) then
         ierr = 1
         err = d1mach(2)
         xernam = name
         write (xern1, '(i8)') locw
         write (xern2, '(i8)') lenw
         call xermsg ('slatec', 'dchkw',
     $      'in ' // xernam // ', double precision work array too ' //
     $      'short.  rwork needs ' // xern1 // '; have allocated ' //
     $      xern2, 1, 1)
      endif
      return
c------------- last line of dchkw follows ----------------------------
      end
