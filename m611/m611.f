*deck @(#)m611.f	5.1 11/6/94
      program m611
c***begin prologue     m611.f
c***date written       940204   (yymmdd)  
c***revision date      11/6/94      
c
c***keywords           dft, density-functional 
c***author             martin, r.l., russo,t.v., hay, p.j. (lanl) 
c***source             @(#)m611.f	5.1   11/6/94
c***purpose            evaluates properties associated with the density
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       m611.f
      implicit none
c     --- local variables ---
      integer ia,ioff,maxcor,intoff,wpadti
      real*8 z(1)
c
      common // ia(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(0,z,maxcor,'m611',0)
      intoff=wpadti(ioff)
      call pm611(z(ioff),ia(intoff))
c
c
      stop
      end
