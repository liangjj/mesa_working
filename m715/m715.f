*deck @(#)m715.f	5.1  11/28/95
      program m715
c***begin prologue     m715.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin, richard and saxe, paul(lanl)
c***source             @(#)m715.f	5.1   11/28/95
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       m715.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer ia,ioff,intoff,maxcor,wpadti
      real*8 z(1)
c
      common // ia(1)
      common /memory/ ioff
      equivalence (ia(1),z(1))
      integer mynodeid, nodeid

c
c
      call pbeginf
c      call niceftn(10)
      if (nodeid() .eq. 0) call drum
c      call  llog
      call getscm(0,z,maxcor,'m715',0)
      intoff=wpadti(ioff)
      call pm715(z(ioff),ia(intoff))
      call pend
c
c
      stop
      end
