*deck @(#)m511.f	5.1 11/6/94
      program m511
c***begin prologue     m511.f
c***date written       930515   (yymmdd)  
c***revision date      11/6/94      
c
c***keywords           kohn-sham, dft, density-functional 
c***author             martin, r.l., russo,t.v., hay, p.j. (lanl) 
c***source             @(#)m511.f	5.1   11/6/94
c***purpose            solves the molecular kohn-sham equations 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       m511.f
      implicit none
c     --- local variables ---
      integer ia,ioff,need,maxcor,intoff,wpadti
      real*8 z(1)
c
      common // ia(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'m511',0)
      intoff=wpadti(ioff)
      call pm511(z(ioff),ia(intoff))
c
c
      stop
      end
