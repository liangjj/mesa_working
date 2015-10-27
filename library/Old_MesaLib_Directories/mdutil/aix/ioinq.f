*deck @(#)ioinq.f	5.1  11/6/94
      function ioinq(file,query)
c
c***begin prologue     ioinq
c***date written       870706   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           iosys dependent routines
c***author             saxe, paul (lanl)
c***source             @(#)ioinq.f	5.1   11/6/94
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
c
c
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
