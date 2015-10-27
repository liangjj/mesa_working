*deck @(#)m250.f	5.1  11/6/94
      program m250
c***begin prologue     m250.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)m250.f	5.1   11/6/94
c***purpose            follows the reaction path hamiltonian
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       m250.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer ia,ioff,intoff,wpadti,maxcor
      real*8 z(1)
c
      common // ia(1)
      common /memory/ ioff
      equivalence (ia(1),z(1))
c
c
      call drum
      call getscm(0,z,maxcor,'m250',0)
      intoff=wpadti(ioff)
      call izero(ia(intoff),maxcor-1)
      call pm250(z(ioff),ia(intoff))
c
c
      stop
      end
