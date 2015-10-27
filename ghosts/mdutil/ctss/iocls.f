*deck %W%  %G%
      integer function iocls(unit,status)
c
c***begin prologue     iocls
c***date written       870706   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           iosys dependent routine
c***author             saxe, paul (lanl)
c***source             %W%   %G%
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
c
      character*(*) status
      integer unit
c
      if (status.eq.'delete') then
         call drabsf(unit,0,-1,ierr)
      else
         call close(unit)
         iocls=0
      end if
c
c
      return
      end
