*deck @(#)m1991.f	5.1  11/6/94
      program m1991
c***begin prologue     m1991.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           population analysis, properties
c***author             martin, richard and saxe, paul(lanl)
c***source             @(#)m1991.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       m1991.f
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
      call drum
      call getscm(0,z,maxcor,'m1991',0)
      intoff=wpadti(ioff)
      call pm1991(z(ioff),ia(intoff))
c
      stop
      end
