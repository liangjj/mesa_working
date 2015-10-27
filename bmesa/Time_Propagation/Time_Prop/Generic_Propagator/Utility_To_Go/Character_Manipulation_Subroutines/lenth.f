*deck  @(#)length.f	1.1 8/2/91
      function lenth(string)
c
c***begin prologue     lenth
c***date written       870715   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)length.f	1.1   11/30/90
c
c***purpose            find the length of a string, without trailing
c                      blanks.
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue     lenth
c
      implicit integer(a-z)
      integer lenth
c
      character*(*) string
c
      do 1 i=len(string),1,-1
         if (string(i:i).ne.' ') then
            lenth=i
            return
         end if
 1    continue
c
c
      lenth=0
      return
      end
