*deck @(#)ioinq.f	5.1  11/6/94
      function ioinq(file,query)
c
c***begin prologue     ioinq
c***date written       870706   (yymmdd)
c***revision date      921007   (yymmdd)
c    7 october 1992    rlm at lanl
c      changing the inquire call to return a variable 'answer'
c      instead of the function ioinq itself.  this apparently caused
c      a compiler bug, when the -static option was invoked.
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
      logical answer
c
      character*(*) file
      character*(*) query
      common/io/inp,iout
c
c
c
      answer=.false.
      if (query.eq.'open') then
         inquire (file=file,opened=answer)
      else if (query.eq.'exist') then
         inquire (file=file,exist=answer)
      end if
      ioinq=answer
c
c
      return
      end
