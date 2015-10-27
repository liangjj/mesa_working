*deck %W%  %G%
      logical function ioinq(file,query)
c
c***begin prologue     ioinq
c***date written       870706   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           iosys dependent routines
c***author             saxe, paul (lanl)
c***source             %W%   %G%
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
c
      character*(*) file
      character*(*) query
      logical temp
c
c
      if (query.eq.'open') then
         temp=.false.
         inquire (file=file,opened=temp)
         ioinq=temp
      else if (query.eq.'exist') then
         temp=.false.
         inquire (file=file,exist=temp)
         ioinq=temp
      end if
c
c
      return
      end
