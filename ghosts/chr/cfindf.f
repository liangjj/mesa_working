*deck @(#)cfindf.f	1.1  11/30/90
      function cfindf(string,substr)
c***begin prologue     cfindf
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           character, string, substring
c***author             martin, richard (lanl)
c***source             @(#)cfindf.f	1.1   11/30/90
c***purpose            finds the first occurrence of a pattern in a string.
c***description
c                      cfindf is an integer function used as:
c                        pos=cfindf(string,substr)
c                          string    the string to search.
c                          substr    the pattern for which to search.
c
c                      cfindf returns the index of the beginning of the
c                      first occurrence of the pattern. it returns
c                      0 if the pattern is not found.
c
c***references
c***routines called    (none)
c***end prologue       cfindf
      implicit integer(a-z)
      integer cfindf
      character*(*) string,substr
c
c
      lenstr=len(string)
      lensub=len(substr)
      cfindf=0
      if(lensub.gt.lenstr) return
      cfindf=index(string,substr)
c
c
      return
      end
