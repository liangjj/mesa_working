*deck %W% %G%
      program xm999
c***begin prologue     %M%
c***date written       940104   (yymmdd)  
c***revision date      %G%
c
c***keywords           least squares
c***author             russo,t.v.
c***source             %M%
c***purpose            solves the least squares problem to fit polynomial to a set of data: version 2, 2D fit
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
      call getscm(need,z,maxcor,'xm999',0)
      intoff=wpadti(ioff)
      call pm999(z(ioff),ia(intoff))
c
c
      stop
      end
