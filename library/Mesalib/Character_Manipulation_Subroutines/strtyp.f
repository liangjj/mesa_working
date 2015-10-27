*deck @(#)strtyp.f	5.1  11/6/94
      function strtyp(string)
c
c***begin prologue     strtyp
c***date written       860815  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           character string type
c***author             saxe, paul (lanl)
c***source             @(#)strtyp.f	5.1   11/6/94
c***purpose            determine whether a string represents an integer,
c                      real*8 number or other character string.
c***description
c               character=strtyp(string)
c
c   on entry:
c             string     character*(*)
c                        the string to be typed.
c
c   on return:
c             strtyp     character*(*)
c                        'integer', 'real', or 'string'
c
c***references
c***routines called    none
c***end prologue       strtyp
c
      implicit integer (a-z)
c
      character*(*) strtyp
      character*(*) string
      character*13 numbs
      character*6  expon
c
c     ----- note the inclusion of blanks in these two strings
c
      parameter (numbs=' +-0123456789',expon=' .eedd')
c
c     ----- check for type of object -----
c
      strtyp='string'
      do 1 i=1,len(string)
         if (index(numbs,string(i:i)).eq.0) go to 2
    1 continue
      strtyp='integer'
      return
c
    2 continue
c
c      ----- if starts with d or e consider as string -----
c
       if (index('eedd',string(1:1)).ne.0) return
c
      do 3 i=1,len(string)
         if (index(numbs//expon,string(i:i)).eq.0) return
    3 continue
      strtyp='real'
c
c
      return
      end
