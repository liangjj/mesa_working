*deck @(#)rotz.f	5.1 11/6/94
      subroutine rotz(t,alpha)
c***begin prologue     rotz.f
c***date written       930711  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin, r.l. (lanl)
c***source             @(#)rotz.f	5.1   11/6/94
c***purpose            returns in t the rotation matrix which rotates the
c                      coordinates counterclockwise about the z-axis by an
c                      angle alpha.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       rotz.f
      implicit none
c     --- input variables -----
      real*8 t(3,3)
      real*8 alpha
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      real*8 one
      parameter (one=1.0d+00)
c
c
      call rzero(t,9)
      t(3,3)=one
      t(1,1)= cos(alpha)
      t(1,2)=-sin(alpha)
      t(2,1)= sin(alpha)
      t(2,2)= cos(alpha) 
c
c
      return
      end
