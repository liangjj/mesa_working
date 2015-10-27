*deck mkh1d.f
c***begin prologue     mkh1d
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            right hand side of inhomogeneous time-dependent
c***                   hamiltonian for 1-dimensional spatial hamiltonian.
c***                   
c***references         
c
c***routines called    
c***end prologue       mkh1d
      subroutine mkh1d(hx,vecout,vecin,tmp,nx)
      implicit integer (a-z)
      real*8 hx, vecout, vecin, tmp
      character*80 title
      dimension hx(nx,nx), vecin(nx), vecout(nx)
      dimension tmp(nx)
      common/io/inp, iout
      call ebc(tmp,hx,vecin,nx,nx,1)
      do 10 j=1,nx
         vecout(j) = vecout(j) + tmp(j)
 10   continue   
      return
      end       
