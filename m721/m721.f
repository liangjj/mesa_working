*deck @(#)m721.f	5.2 2/5/95
      program m721
c***begin prologue     m721.f
c***date written       yymmdd  
c***revision date      2/5/95      
c
c***keywords           
c***author             
c***source             @(#)m721.f	5.2   2/5/95
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       m721.f
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
      call getscm(0,z,maxcor,'m721',0)
      intoff=wpadti(ioff)
      call pm721(z(ioff),ia(intoff))
c
c
      stop
      end
