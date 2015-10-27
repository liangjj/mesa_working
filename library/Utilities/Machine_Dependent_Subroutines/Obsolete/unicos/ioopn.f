*deck @(#)ioopn.f	5.1  11/6/94
      subroutine ioopn(unit,file,status,iostat,extra)
c
c***begin prologue     ioopn
c***date written       870706   (yymmdd)
c***revision date      910617   (yymmdd)
c
c   71 june   1991     rlm at lanl
c   unicos version
c***keywords           iosys dependent routine
c***author             saxe, paul (lanl)
c***source             @(#)ioopn.f	5.1   11/6/94
c
c***purpose            to open a unit to a disc file 
c
c***description        
c
c***references         
c
c***routines called    (none)
c
c***end prologue       ioopn
c
      implicit integer (a-z)
c
      parameter (maxun=20)
c
      character*(*) file
      character*(*) status
      character*(*) extra
      character*4 itoc
      character*8 iounq
      character*256 unitnm
      integer unit
      integer iostat
      common/io/inp,iout
      common /fnames/ unitnm(maxun)
c
c     ----- open the unit as a fortran unit -----
c
      if (status.eq.'scratch') then
c        generate a unique filename for scratch files.
c        pass this back to ioopen.
         file=iounq()
         unitnm(unit)=file(1:7)
         call wopen(file(1:7),20,0,ierr)
         if(ierr.ne.0) then
            call lnkerr('problem in wopen:'//itoc(ierr))
         endif
      else
         unitnm(unit)=file(1:7)
         call wopen(file(1:7),20,0,ierr)
         if(ierr.ne.0) then
            call lnkerr('problem in wopen:'//itoc(ierr))
         endif
      end if
c
c
      return
      end
