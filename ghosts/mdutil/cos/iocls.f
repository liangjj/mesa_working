*deck %W%  %G%
      subroutine iocls(unit,status)
c
      implicit integer (a-z)
c
      character*(*) status
c
      call wclose(unit)
c
      if (status.eq.'delete') then
         close (unit=unit,status='delete')
      else
         close (unit=unit)
      end if
c
c
      return
      end
