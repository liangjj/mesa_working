*deck @(#)m6203.f	1.2  10/27/94
      program m6203
c***begin prologue     m6203.f
c***date written       930515   (yymmdd)  
c***revision date      10/27/94      
c
c***keywords           poisson
c***author             schneider, barry(nsf)
c***source             @(#)m6203.f	1.2   10/27/94
c***purpose            solves the generalized molecular poisson equation
c                      del**2 +k**2= rho 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       m6203.f
      implicit none
c     --- local variables ---
      integer ia,ioff,need,maxcor,intoff,wpadti
      real*8 z(1)
c
      common // ia(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'m6203',0)
      intoff=wpadti(ioff)
      call pm6203(z(ioff),ia(intoff))
c
c
      stop
      end
