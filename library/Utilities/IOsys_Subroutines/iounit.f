*deck @(#)iounit.f	5.1  11/6/94
      function iounit(unit,unlist,nunits,error)
c
c***begin prologue     iounit
c***date written       850125   (yymmdd)
c***revision date      910225   (yymmdd)
c
c    february 25,1991  rlm at lanl
c       renaming loop counter so that it doesn't conflict with the function
c       name -- a problem on the alliant.
c
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)iounit.f	5.1   11/6/94
c***purpose            to return the unit number of an iosys unit.
c
c***description        #
c
c
c***references
c
c***routines called    (none)
c
c   common blocks:     (none)
c
c***end prologue       iounit
c
      implicit integer (a-z)
      integer iounit
c
      character*(*) unit,unlist(nunits)
c
      iounit=0
      error=0
      do 1 i=1,nunits
         if (unit.eq.unlist(i)) then
            iounit=i
            return
         end if
    1 continue
c
c     ----- error condition exists -----
c
      error=1
c
      end
