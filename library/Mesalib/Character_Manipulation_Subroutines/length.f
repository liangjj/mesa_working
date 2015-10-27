*deck  @(#)length.f	1.1 8/2/91
      function length(string)
c
c***begin prologue     length
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
c***end prologue     length
c
      implicit integer(a-z)
      integer length, lenth
      character*(*) string
      length=lenth(string)
c
      return
      end
