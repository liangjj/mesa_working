*deck %W%  %G%
      subroutine ioopn(unit,file,status,iostat,extra)
c
c***begin prologue     ioopn
c***date written       870706   (yymmdd)
c***revision date      891219   (yymmdd)
c
c                      rlm at lanl
c                      the size parameter which is passed through extra
c                      is now an estimate of the total file size, as 
c                      opposed to the family member length.  we
c                      have a maximum of 99 family members at lanl, and
c                      this limit is reflected in the call to famsiz below.
c                      the minimum family member length is set to 4mwords.
c 
c***keywords           iosys dependent routine
c***author             saxe, paul (lanl)
c***source             %W%   %G%
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
      character*(*) file
      character*(*) status
      character*(*) extra
      character*8 iounq
      character*8 scrnam
      integer unit
      integer intkey
      logical logkey
c
c     ----- open the unit as a cftlib random, familied unit -----
c
      if (status.eq.'scratch') then
c
c        ----- get a unique file-name -----
c
         scrnam=iounq()
         call assign(unit,scrnam,0)
      else
         call assign(unit,file,0)
      end if
c
c     ----- put the file on the ssd if requested -----
c
      if (logkey(extra,'ssd',.false.,' ')) then
         call setssd(unit)
      end if
c
c     ----- set the family size and wait parameters -----
c
      size=intkey(extra,'size',262144,' ')
      size=max(size/99,4000000)
      call famsiz(unit,size)
      call famwait(unit,1)
c
c
      return
      end
