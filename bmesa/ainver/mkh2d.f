*deck mkh2d.f
c***begin prologue     mkh2d
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            right hand side of inhomogeneous time-dependent
c***                   hamiltonian for 2-dimensional spatial hamiltonian.
c***                   
c***references         
c
c***routines called    
c***end prologue       mkh2d
      subroutine mkh2d(hx,hy,driver,psi0,tmp,nx,ny)
      implicit integer (a-z)
      real*8 hx, hy, vecout, vecin, tmp
      dimension hx(nx,nx), hy(ny,ny), vecout(ny,nx)
      dimension vecin(ny,nx), tmp(ny,nx)
      common/io/inp, iout
      call ebc(tmp,hy,vecin,ny,ny,nx)
      do 10 j=1,nx
         do 20 k=1,ny
            vecout(k,j) = vecout(k,j) + tmp(k,j) 
 20      continue   
 10   continue
      call ebc(tmp,vecin,hx,ny,nx,nx)
      do 100 j=1,nx
         do 200 k=1,ny
            vecout(k,j) = vecout(k,j) + tmp(k,j) 
 200     continue   
 100  continue   
      return
      end       
