*deck @(#)laplace.f	5.1 11/6/94
      subroutine laplace(r,j,y,nr,lmax)
c***begin prologue     laplace.f
c***date written       940304    (yymmdd)  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin, richard(lanl) and schneider, barry(nsf)
c***source             @(#)laplace.f	5.1   11/6/94
c***purpose            
c***description
c   this sets up the homogeneous solutions: 
c      r**(l+1) in j, and -(1/r)**l in y for all l from 0:lmax.
c
c      note that for the poisson equation, these are the green's functions 
c      for (r*rho), which yields a function u=r*v
c
c***references
c
c***routines called
c
c***end prologue       laplace.f
      implicit none
c     --- input variables -----
      integer nr,lmax
c     --- input arrays (unmodified) ---
      real*8 r(nr)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 j(nr,0:lmax),y(nr,0:lmax)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer pt,l
      real*8 zero,one
c
      parameter (zero=0.0d+00,one=1.0d+00)
      common/io/ inp,iout
c
      do 20 pt=1,nr
         do 10 l=0,lmax
            j(pt,l)=r(pt)**(l+1)
            if(pt.eq.1) then
               y(pt,l)=zero
            else
               y(pt,l)=-(one/r(pt))**l
            endif
   10    continue
   20 continue
c
c
      return
      end
