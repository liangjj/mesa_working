*deck @(#)iocls.f	4.1  7/7/93
      function iocls(unit,status)
c
c***begin prologue     iocls
c***date written       870706   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           iosys dependent routine
c***author             saxe, paul (lanl)
c***source             @(#)iocls.f	4.1   7/7/93
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
      character*(*) status
      integer unit
c
      close (unit=unit,status=status,iostat=iocls)
c
c
      return
      end
