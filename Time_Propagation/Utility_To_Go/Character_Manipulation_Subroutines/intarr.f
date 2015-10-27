*deck @(#)intarr.f	5.1  11/6/94
      subroutine intarr(string,key,array,n,ref)
c
c***begin prologue     intarr
c***date written       870901  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c
c***keywords           integer keywords, namelist
c***author             saxe, paul (lanl)
c***source             @(#)intarr.f	5.1   11/6/94
c***purpose            to return an array of integers associated
c     with a keyword.
c***description
c               call intarr(string,key,array,n,ref)
c
c    intarr searches string for a nested set of keywords in key and returns
c    either the integer value associated with  the keyword if the keyword
c    is found.
c
c   on input:
c             string      character*(*)
c                         the string to search for the keyword(s).
c
c             key         character*(*)
c                         the keyword(s) to search for.
c
c             array       integer(n)
c                         the array of values to fill in.
c
c             n           integer.
c                         the maximum number of elements.
c
c             ref         character*(*)
c                         the reference string of keywords for truncation
c                         searches. a blank string disables matching.
c
c
c
c***references
c***routines called    chrkey      (char)
c                      ctoi        (char)
c***end prologue       intarr
c
      implicit integer (a-z)
c
      integer array(n)
      character*(*) string,key,ref
      character*20 integr
      character*1 nxtkey
c
c     -----
c
      call poskey(string,key,ref,start,end)
c
c     ----- blank indicates couldn't find the string -----
c
      if (start.le.0) then
         return
      end if
c
c     ----- check for non-integer string, and return default for all such
c           cases
c
      ndone=0
      pos=0
 1    continue
         if (nxtkey(string(start:end),pos,integr).eq.' ') go to 10
         ndone=ndone+1
         if (ndone.gt.n) then
            integr=key
            call lnkerr('too many values supplied to intarr '//integr)
         end if
c
         array(ndone)=ctoi(integr)
      go to 1
c
 10   continue
c
c
      return
      end
