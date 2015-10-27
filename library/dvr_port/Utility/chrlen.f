*deck  @(#)chrlen.f	1.1 8/2/91
      function chrlen(string)
c
c***begin prologue     chrlen
c***date written       870715   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)chrlen.f	1.1   11/30/90
c
c***purpose            find the length of a string, without trailing
c                      blanks.
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue     chrlen
c
      implicit integer(a-z)
      integer chrlen
c
      character*(*) string
c
      do 1 i=len(string),1,-1
         if (string(i:i).ne.' ') then
            chrlen=i
            return
         end if
 1    continue
c
c
      chrlen=0
      return
      end
