*deck @(#)lnkerr.f	5.1  11/6/94
      subroutine lnkerr(messge)
c***begin prologue     lnkerr
c***date written       850601  (yymmdd)
c***revision date      870207  (yymmdd)
c
c   7 february 1987   pws at lanl
c         modifying for bsd 4.2 unix on sun 3/50 and 3/160 workstations
c
c   1 december 1986   pws at lanl
c         modifying to close the input and output files
c
c***author             martin, richard (lanl)
c***source             @(#)lnkerr.f	5.1   11/6/94
c***purpose            error termination routine for mesa.
c***description
c                      this routine is used for error exits from mesa.
c                      it first writes the message to the print file,
c                      generates subroutine traceback information, and
c                      aborts the job.
c
c                      call lnkerr(messge)
c                        messge  an informative message which is printed
c                                before shutdown. only the first 123
c                                characters will be printed. character*(*).
c                                if the message contains the string 'noclose',
c                                then the iosys files are not closed before
c                                aborting.  this is used for lnkerr calls from
c                                within iosys.  if the error occurs in iosys,
c                                an attempt to close the files may generate
c                                recursive calls.
c
c***references
c***routines called    ioabor(io), trcbk(ctss), abort(ctss)
c***end prologue       lnkerr
c
      implicit integer(a-z)
      character*(*) messge
      common/io/inp,iout
c
 1000 format(' lnkerr: ',123a1)
c
      write(iout,1000) (messge(i:i),i=1,len(messge))
      if(index(messge,'noclose').eq.0) then
         call ioabor
      endif
      call trbk(iout,0)
c
c     close the external communication files.
c
      close (inp)
      close (iout)
c
      call abort(messge)
c
c
      stop
      end
