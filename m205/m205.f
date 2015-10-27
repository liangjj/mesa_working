*deck @(#)m205.f	5.1  11/6/94
      program m205
c***begin prologue     m205.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)m205.f	5.1   11/6/94
c***purpose            computes vibrational frequencies from 
c                      finite difference second derivatives.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       m205.f
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
      call getscm(0,z,maxcor,'m205',0)
      intoff=wpadti(ioff)
c     --- initialize core to zero.
      call izero(ia(intoff),maxcor-1)
      call pm205(z(ioff),ia(intoff))
c
c
      stop
      end
