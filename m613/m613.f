*deck @(#)m613.f	5.1  11/6/94
      program m613
c***begin prologue     m613.f
c***date written       940304   (yymmdd)  
c***revision date      11/6/94      
c
c***keywords           poisson
c***author             martin, richard(lanl) and schneider, barry(nsf)
c***source             @(#)m613.f	5.1   11/6/94
c***purpose            solves the generalized molecular poisson equation
c                      (del**2 +k**2)v= rho 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       m613.f
      implicit none
c     --- local variables ---
      integer ia,ioff,maxcor,intoff,wpadti,need
      real*8 z(1)
c
      common // ia(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'m613',0)
      intoff=wpadti(ioff)
      call pm613(z(ioff),ia(intoff))
c
c
      stop
      end
