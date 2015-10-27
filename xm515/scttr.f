*deck @(#)scttr.f	1.1  4/25/95
      subroutine scttr(rout,rin,inds,n)
c***begin prologue     scttr.f
c***date written       930618  
c***revision date      11/6/94      
c
c***keywords           scatter
c***author             RUSSO, thomas (lanl)
c***source             @(#)scttr.f	1.1   4/25/95
c***purpose            scatter rin into rout using inds
c                       
c                      
c                      
c***description
c
c     real*8 rin(*)      input array to scatter from
c     real*8 rout(n)      output array to scatter into
c     integer i(n)     array of indices
c     integer n        number of things to gather
c
c***references         
c none
c
c***routines called
c
c***end prologue       scttr.f
      implicit none
c --input variables-- (unmodified)      
      integer n
c --input arrays-- (unmodified)
      real*8 rin(n)
      integer inds(n)
c --output variables--
      real*8 rout(*)
c local
      integer i
      do 10 i=1,n
         rout(inds(i))=rin(i)
 10   continue 
      return
      end
