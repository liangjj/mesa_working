*deck @(#)cskipf.f	5.1  11/6/94
      function cskipf(string,substr)
c***begin prologue     cskipf
c***date written       850601  (yymmdd)
c***revision date      910304  (yymmdd)
c
c    3 march, 1991     rlm at lanl
c      defining the function early on. there was the possibility of
c      a return before it was defined. 
c***keywords           character, string, skip
c***author             martin, richard (lanl)
c***source             @(#)cskipf.f	5.1   11/6/94
c***purpose            skips forward over all occurrences of a pattern.
c***description
c                      cskipf is an integer function used as:
c                        pos=cskipf(string,substr)
c                          string  input character string to search.
c                          substr  substring to skip over.
c
c                      cskipf returns the index of the first character
c                      which does not match the pattern.  it returns 0
c                      if the string contains only the pattern of interest.
c
c***references
c***routines called    (none)
c***end prologue       cskipf
      implicit integer(a-z)
      integer cskipf
      character*(*) string,substr
c
c
c
      cskipf=0
      lenstr=len(string)
      lensub=len(substr)
      if(lensub.gt.lenstr) return
      do 10 i=1,lenstr,lensub
         if(string(i:i-1+lensub).ne.substr) then
            cskipf=i
            return
         endif
   10 continue
c
c
      return
      end
