*deck @(#)nxtlin.f	5.1  11/6/94
      function nxtlin(line,unit,lineno)
c***begin prologue     nxtlin
c***date written       850601  (yymmdd)
c***revision date      870207  (yymmdd)
c
c   7 february 1987   pws at lanl
c     making a character function.
c
c***keywords           input/output
c***author             saxe, paul (lanl).
c***source
c***purpose            reads the next line in an input file.
c***description
c                      nxtlin is a character function used as:
c                      line=nxtlin(line,unit,lineno)
c
c                        line      the next line (character*80)
c                        unit      the unit number of the device to search.
c                        lineno    a variable which is incremented for each
c                                  line read.
c
c                      nxtlin='ok'  if a line was successfully read, and
c                      nxtlin='eof' if an end-of-file was encountered.
c
c                      lines beginning with blanks or '/' are skipped,
c                      and the text is capitalized before the return.
c
c***references
c***routines called    locase(chr)
c***end prologue       nxtlin
      implicit integer(a-z)
      character*(*) nxtlin
      character*80 line
      integer unit
c
    1 continue
         read (unit,2,end=9) line
    2    format (a80)
         lineno=lineno+1
      if (line.eq.' '.or.line(1:1).eq.'/') go to 1
c
      nxtlin='ok'
      call locase(line,line)
      return
c
    9 continue
      nxtlin='eof'
      return
      end
