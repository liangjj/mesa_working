*deck @(#)str.f	5.1  11/6/94
      subroutine str(noint,i,j,b,ib,c)
c***begin prologue     str.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)str.f	5.1   11/6/94
c***purpose            
c***description
c     
c        adapted from the normal coordinate analysis program of
c        schachtschneider, shell development .
c
c***references
c
c***routines called
c
c***end prologue       str.f
      implicit none
c     --- input variables -----
      integer noint,i,j
c     --- input arrays (unmodified) ---
      real*8 c(*)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer ib(4,2)
      real*8 b(3,4,2)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer iaind,jaind,m
      real*8 rij(3)
      real*8 zero,dijsq
      parameter (zero=0.0d+00)
c
c
      iaind=3*(i-1)
      jaind=3*(j-1)
      ib(1,noint)=i
      ib(2,noint)=j
      dijsq = zero
      do 10 m=1,3
         rij(m)=c(m+jaind)-c(m+iaind)
         dijsq=dijsq+rij(m)**2
   10 continue
      do 20 m=1,3
         b(m,1,noint)=-rij(m)/sqrt(dijsq)
         b(m,2,noint)=-b(m,1,noint)
   20 continue
c
c
      return
      end
