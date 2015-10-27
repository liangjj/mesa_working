*deck @(#)gthr.f	5.2 2/5/95
      subroutine gthr(rout,rin,inds,n)
c***begin prologue     gthr.f
c***date written       930618  
c***revision date      2/5/95      
c
c***keywords           gather
c***author             RUSSO, thomas (lanl)
c***source             @(#)gthr.f	5.2   2/5/95
c***purpose            gather stuff from rin into rout using inds
c                       
c                      
c                      
c***description
c
c     real*8 rin(*)    array from which to gather nonzeros
c     real*8 rout(n)   array into which to gather nonzeros
c     integer i(n)     array of indices
c     integer n        number of things to gather
c
c***references         
c
c***routines called
c
c***end prologue       gthr.f
      implicit none
c     --- input variables -----
      integer n
c     --- input arrays (unmodified) ---
      real*8 rin(*)
      integer inds(n)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 rout(n)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i
c
c
      do 10 i=1,n
         rout(i)=rin(inds(i))
 10   continue 
c
c
      return
      end
