*deck @(#)cskipb.f	5.1  11/6/94
      function cskipb(string,substr)
c***begin prologue     cskipb
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           character, string, skip
c***author             martin, richard (lanl)
c***source             @(#)cskipb.f	5.1   11/6/94
c***purpose            skips backwards over all occurrences of a pattern.
c***description
c                      cskipb is an integer function used as:
c                        pos=cskipb(string,substr)
c                          string  input character string to search.
c                          substr  string to skip over.
c
c                      cskipb returns the index of the last character
c                      which does not match the pattern.  it returns 0
c                      if the string contains only the pattern of interest.
c                      it is generally used to return the length of a string
c                      ignoring trailing blanks, i.e; length=cksipb(str,' ')
c***references
c***routines called    (none)
c***end prologue       cskipb
      implicit integer(a-z)
      integer cskipb
      character*(*) string,substr
c
c
      lenstr=len(string)
      lensub=len(substr)
      cskipb=0
      if(lensub.gt.lenstr) return
c
      do 10 i=lenstr,1,-lensub
         if(string(i+1-lensub:i).ne.substr) then
            cskipb=i
            return
         endif
   10 continue
c
c
      return
      end
