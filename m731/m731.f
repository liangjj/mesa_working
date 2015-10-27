*deck @(#)m731.f	5.1  11/6/94
      program m731
c***begin prologue     m731.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)m731.f	5.1   11/6/94
c***purpose            transforms forces from cartesians to internals.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       m731.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer ia
      integer maxcor,ioff,intoff,wpadti
      real*8 z(1)
c
      common // ia(1)
      common /memory/ ioff
      equivalence (ia(1),z(1))
c
c
      call drum
      call getscm(0,z,maxcor,'m731',0)
      intoff=wpadti(ioff)
      call pm731(z(ioff),ia(intoff))
c
c
      stop
      end
