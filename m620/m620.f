*deck @(#)m620.f	5.1  11/6/94
      program m620
c***begin prologue     m620.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin, richard (lanl)
c***source             @(#)m620.f	5.1   11/6/94
c***purpose            cox-williams potential fit
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       m620.f
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
c
c
      call drum
      call getscm(0,z,maxcor,'m620',0)
      intoff=wpadti(ioff)
      call pm620(z(ioff),ia(intoff))
      call chainx(0)
c
c
      stop
      end
