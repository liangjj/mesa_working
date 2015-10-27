*deck hamfl1.f
c***begin prologue     hamfl1
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form all or a sub-block of the two dimensional 
c***                   hamiltonian.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       hamfl1
      subroutine hamfl1(h,hx,vx,nx,prn)
      implicit integer (a-z)
      real*8 h, hx, vx
      character*80 title 
      logical prn
      dimension h(nx,nx), hx(nx,nx), vx(nx)
      common/io/inp, iout
      do 10 xi=1,nx
         do 20 xj=1,nx
            h(xi,xj) = h(xi,xj) + hx(xi,xj)
 20      continue
 10   continue   
      do 30 xi=1,nx
         h(xi,xi) = h(xi,xi) + vx(xi)
 30   continue   
      return
      end       

