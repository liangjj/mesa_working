*deck @(#)nzident.f	5.1 11/6/94
      subroutine nzident(r,i,n)
c***begin prologue     nzident.f
c***date written       930618  
c***revision date      11/6/94      
c
c***keywords           non-zero indices prepare for gather
c***author             RUSSO, thomas (lanl)
c***source             @(#)nzident.f	5.1   11/6/94
c***purpose            Fills array i with indices of nonzero elements
c                      of array r.  This is to prepare r for a gather 
c                      operation
c                      
c***description
c
c     real*8 r(n)      Array from which to gather nonzeros
c     integer i(n)     array of nonzeros
c     integer n        size of r.  On return will have number of nonzeros
c
c***references         
c  References?  They didn't give us no references.  We don't gotta show you
c  no references.  We don't NEED no stinkin' references.
c
c***routines called
c
c***end prologue       nzident.f

      implicit none

c --input variables--(modified)
      integer n
c --input arrays--    (unmodified)
      real*8 r(n)
c --output arrays--
      integer i(n)
c --local vars--
      integer j,ntmp
      real*8 zero
      parameter (zero=0.0d0)
      ntmp=0
      do 10 j=1,n
         if (r(j) .ne. zero) then
            ntmp=ntmp+1
            i(ntmp)=j
         endif
 10   continue 
      n=ntmp
      return
      end
