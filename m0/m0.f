*deck @(#)m0.f	5.1 11/6/94
      program m0
c***begin prologue     m0.f
c***date written       850601   (yymmdd)
c***revision date      11/6/94
c
c    1 february 1993   rlm at lanl
c        adding a machine-dependent routine to return the unit number
c        for standard error.  this should make this routine less
c        machine-dependent.
c   22 may 1992        rlm at lanl
c        updating the commentary
c    7 february 1987   pws at lanl
c        making modifications to create bsd 4.2 unix version on
c        sun 3/52 and 3/160 workstations.
c
c***keywords           m0, link 0, driver, i/o
c***author             martin, richard  (lanl)
c***source             @(#)m0.f	5.1 11/6/94
c***purpose            main driver for mesa system.
c***description
c     this link is (in principle) the only link is the main sequence
c     which is machine dependent. the call to abort may be machine dependent.
c     at present, this program simply prepares the output file
c     for later use. note that the unit number assignments 
c     of input=8 output=9 must be compatible with the definitions in
c     lxopen(mdutil).
c     the sequence of events is:
c        check for file name replacement
c        destroy any existing file with the same name, and assign
c        a new one.
c        print out the current version of the code
c        close the output fileand notify the controller
c        to start up link1.
c
c
c***references
c***routines called    ioinq(mdutil), iorm(mdutil), versn(mdutil)
c                      comand(mdutil), chain(util)
c***end prologue       m0.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      logical ioinq
      character*128 outnam
      character*2 next
      integer stderr
      integer inp,iout,ierr
c
      common/io/inp,iout
c
c     --- assign input and output unit numbers, and check the command
c         line to see if the output file name was changed from default.
      inp=8
      iout=9
      ierr=stderr()
      outnam='out'
      call comand(1,outnam)
      if (outnam.eq.'out') outnam='mesa.out'
c
c     --- destroy this file if it exists ---
      if (ioinq(outnam,'exist')) then
         call iorm(outnam)
      end if
c
c     --- open the output file ---
      open (unit=iout,file=outnam,access='sequential',form='formatted',
     $      err=1,status='unknown')
c
c     --- write current version of mesa to the output file ---
      call versn(iout)
c
c     --- close the output file and send the next link name to the
c         controller
      close (iout)
      next='m1'
      call chain(next)
c
c
c      stop
      Call exit
c
c     --- handle an inability to open the main output file ---
    1 continue
      write(ierr,*) 'link 0 cannot open the output file'
      call abort()
c
c
      end
