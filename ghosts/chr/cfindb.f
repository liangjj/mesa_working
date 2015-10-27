*deck @(#)cfindb.f	1.1  11/30/90
      function cfindb(string,substr)
c***begin prologue     cfindb
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           character, string, substring
c***author             martin, richard (lanl)
c***source             @(#)cfindb.f	1.1   11/30/90
c***purpose            finds the last occurrence of a pattern in a string.
c***description
c                      cfindb is an integer function used as:
c                        pos=cfindb(string,substr)
c                          string    the string to search.
c                          substr    the pattern for which to search.
c
c                      cfindb returns the index of the beginning of the
c                      last occurrence of the pattern. it returns
c                      0 if the pattern is not found.
c
c***references
c***routines called    (none)
c***end prologue       cfindb
      implicit integer(a-z)
      integer cfindb
      character*(*) string,substr
c
c
      lenstr=len(string)
      lensub=len(substr)
      cfindb=0
      if(lensub.gt.lenstr) return
      do 10 i=lenstr,lensub,-1
         if(string(i+1-lensub:i).eq.substr) then
            cfindb=i+1-lensub
            return
         endif
   10 continue
c
c
      return
      end
