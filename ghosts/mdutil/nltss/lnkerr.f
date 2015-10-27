*deck %W%  %G%
      subroutine lnkerr(messge)
c***begin prologue     lnkerr
c***date written       850601  (yymmdd)
c***revision date      861201  (yymmdd)
c
c   1 december 1986   pws at lanl
c         modifying to close the input and output files
c
c***author             martin, richard (lanl)
c***source             %W%   %G%
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
      implicit integer(a-z)
      character*(*) messge
      common/io/inp,iout
c
 1000 format(' lnkerr: ',123a1)
c
      write(iout,1000) (messge(i:i),i=1,len(messge))
      if(index(messge,'noclose').eq.0) then
c        call trakio(0,-1)
         call ioabor
      endif
c     call trcbk(iout,0)
c
c     close the external communication files.
c
      call close(inp)
      call close(iout)
c
      call abort('lnkerr.')
      stop
      end
