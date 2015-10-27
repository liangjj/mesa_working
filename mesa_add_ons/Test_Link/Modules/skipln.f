*deck @(#)skipln.f	5.1  11/6/94
      function skipln(unit,n)
c
c***begin prologue     skipln
c***date written       850601  (yymmdd)
c***revision date      870207  (yymmdd)
c
c   7 february 1987   pws at lanl
c     making skipln a character function.
c
c***keywords           input/output, skip
c***author             saxe, paul (lanl).
c***source
c***purpose            skips n lines in an input file.
c***description
c                      skipln is an character function used as:
c                      line=skipln(unit,n)
c
c                        unit      the unit number of the device to search.
c                        n         the number of lines to skip.
c
c                      skipln='ok'  if n lines were successfully skipped, and
c                      skipln='eof' if an end-of-file was encountered.
c
c
c***references
c***routines called    (none)
c***end prologue       skipln
      implicit integer(a-z)
      character*(*) skipln
      integer unit,n
c
c
      rewind unit
      do 1 i=1,n
         read (unit,*,end=2)
    1 continue
      skipln='ok'
      return
c
    2 continue
      skipln='eof'
c
c
      return
      end
