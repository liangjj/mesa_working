*deck @(#)m621.f	5.1  11/6/94
      program m621
c***begin prologue     m621.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             tawa, greg and martin, richard(lanl)
c***source             @(#)m621.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       m621.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer ia,wpadti,ioff,intoff,maxcor
      real*8 z(1)
c
      common // ia(1)
      common /memory/ ioff
      equivalence (ia(1),z(1))
c
c
      call drum
      call getscm(0,z,maxcor,'m621',0)
      intoff=wpadti(ioff)
      call pm621(z(ioff),ia(intoff))
c
c
      stop
      end
