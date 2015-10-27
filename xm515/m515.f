*deck @(#)m515.f	1.2  11/28/95
      program m515
c***begin prologue     m515.f
c***date written       930515   (yymmdd)  
c***revision date      11/6/94      
c
c***keywords           kohn-sham, dft, density-functional 
c***author             martin, r.l., russo,t.v., hay, p.j. (lanl) 
c***source             @(#)m515.f	1.2   11/28/95
c***purpose            solves the molecular kohn-sham equations in parallel
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       m515.f
      implicit none
c     --- local variables ---
      integer ia,ioff,need,maxcor,intoff,wpadti
      real*8 z(1)
      integer nodeid
c
      common // ia(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call pbeginf
      call niceftn(10)
      if (nodeid() .eq. 0) call drum
c debugging      call  llog
      call getscm(need,z,maxcor,'m515',0)
      intoff=wpadti(ioff)
      call pm515(z(ioff),ia(intoff))
      call stats
      call pend
c
c
      stop
      end
