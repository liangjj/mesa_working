*deck @(#)iocls.f	5.1  11/6/94
      function iocls(unit,status)
c
c***begin prologue     iocls
c***date written       870706   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           iosys dependent routine
c***author             saxe, paul (lanl)
c***source             @(#)iocls.f	5.1   11/6/94
c
c***purpose            to close a unit.
c
c***description        
c
c***references         
c
c***routines called    (none)
c
c***end prologue       iocls
c
      implicit integer (a-z)
      integer iocls
c
      parameter (maxun=20)
c
      character*(*) status
      integer unit
c
      character*4 itoc
      character*256 unitnm
      character*8 filenm
c
      common /fnames/ unitnm(maxun)
c
      iocls=0
      filenm=unitnm(unit)
      if(status.eq.'delete') then
         call wclose(filenm,ierr)
         if(ierr.ne.0) then
            call lnkerr('error in wclose:'//itoc(ierr))
         endif
         call iorm(filenm)
      else 
         call wclose(filenm,ierr)
         if(ierr.ne.0) then
            call lnkerr('error in wclose:'//itoc(ierr))
         endif
      endif
c
c
      return
      end
