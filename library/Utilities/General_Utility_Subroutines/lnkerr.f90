!deck @(#)linkerr.f90	1.1  9/6/91
      subroutine linkerr(messge)
!***begin prologue     lnkerr
!***date written       850601  (yymmdd)
!***revision date      870207  (yymmdd)
!
!   7 february 1987   pws at lanl
!         modifying for bsd 4.2 unix on sun 3/50 and 3/160 workstations
!
!   1 december 1986   pws at lanl
!         modifying to close the input and output files
!
!***author             martin, ri!hard (lanl)
!***source             @(#)lnkerr.f	1.1   9/6/91
!***purpose            error termination routine for mesa.
!***description
!                      this routine is used for error exits from mesa.
!                      it first writes the message to the print file,
!                      generates subroutine traceback information, and
!                      aborts the job.
!
!                      call lnkerr(messge)
!                        messge  an informative message which is printed
!                                before shutdown. only the first 123
!                                !haracters will be printed. character*(*).
!                                if the message contains the string 'noclose',
!                                then the iosys files are not closed before
!                                aborting.  this is used for lnkerr calls from
!                                within iosys.  if the error occurs in iosys,
!                                an attempt to close the files may generate
!                                recursive calls.
!
!***references
!***routines called    ioabor(io), trcbk(ctss), abort(ctss)
!***end prologue       lnkerr
!
      IMPLICIT NONE
      CHARACTER(LEN=*)       :: messge
      INTEGER                :: inp
      INTEGER                :: iout
      INTEGER                :: i
      INTEGER                :: index
      INTEGER                :: ierr
      INTEGER                :: stderr
      COMMON/io/ inp, iout
!
 1000 Format(' lnkerr: ',123a1)
!
      Write(iout,1000) (messge(i:i),i=1,len(messge))
      if(index(messge,'noclose').eq.0) then
         call ioabor
      endif
!
!     close the external communication files.
!
      Close (inp)
      Close (iout)
      ierr=stderr()
      write(ierr,*) messge
      stop
End Subroutine  linkerr
