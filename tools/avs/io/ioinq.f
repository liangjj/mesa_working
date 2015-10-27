*deck @(#)ioinq.f	4.1  7/7/93
      function ioinq(file,query)
c
c***begin prologue     ioinq
c***date written       870706   (yymmdd)
c***revision date      910625   (yymmdd)
c   25 june   1991     rlm at lanl
c      isolated a bug in the inquire function in version 1.3.1 of f77
c      on the sun.  the 'open' query gives undependable results.
c      Sun is aware of it and it is fixed in version 1.4.1
c      this function is only called with the 'open' query by io/ioopen.f
c      temporarily getting past it by always returning false for the 
c      'open' query. this may cause a problem if the input or output
c      files happen to have the same name as an iosys file, but so it goes.
c      the inquire statement should be reinstated when v1.4.1 of f77 is
c      available.
c
c***keywords           iosys dependent routines
c***author             saxe, paul (lanl)
c***source             @(#)ioinq.f	4.1   7/7/93
c
c***purpose            determine whether a file is opened, or exists.
c
c***description        
c
c***references         
c
c***routines called    (none)
c
c***end prologue       ioinq
c
      implicit integer (a-z)
      logical ioinq
c
      character*(*) file
      character*(*) query
      common/io/inp,iout
c
c
      ioinq=.false.
      if (query.eq.'open') then
         ioinq=.false.
         inquire (file=file,opened=ioinq)
      else if (query.eq.'exist') then
         ioinq=.false.
         inquire (file=file,exist=ioinq)
      end if
c
c
      return
      end
