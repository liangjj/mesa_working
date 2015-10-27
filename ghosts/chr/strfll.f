*deck @(#)strfll.f	1.1  11/30/90
      subroutine strfll(string,substr)
c***begin prologue     strfll
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           string, fill
c***author             martin, richard (lanl)
c***source             @(#)strfll.f	1.1   11/30/90
c***purpose            fills a string with a pattern.
c***description
c                      call strfll(string,substr)
c                        string  the string to fill.
c                        substr  the pattern with which to fill.
c***references
c***routines called    (none)
c***end prologue       strfll
      implicit integer(a-z)
      character*(*) string,substr
c
c
      lenstr=len(string)
      lensub=len(substr)
      if(lenstr.le.0) return
      string=substr
      if(lenstr.gt.lensub) then
         do 10 i=lensub+1,lenstr,lensub
            string(i:)=substr
   10    continue
      endif
c
c
      return
      end
