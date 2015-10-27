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
      common /blname/ bhlnam(100)
      character*8 bhlnam
c
      character*(*) status
      integer unit
c
      if (status.eq.'delete') then
c        call drabs(unit,0,-1,ierr)
      call destroy(bhlnam(unit))
      else
         call close(unit)
         iocls=0
      end if
c
c
      return
      end
