*deck mkh3d.f
c***begin prologue     mkh3d
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            right hand side of inhomogeneous time-dependent
c***                   hamiltonian for 3-dimensional spatial hamiltonian.
c***                   
c***references         
c
c***routines called    
c***end prologue       mkh3d
      subroutine mkh3d(hx,hy,hz,vecout,vecin,tmp,nx,ny,nz)
      implicit integer (a-z)
      real*8 hx, hy, hz, vecout, vecin, tmp
      dimension hx(nx,nx), hy(ny,ny), hz(nz,nz), vecout(nz,ny,nx)
      dimension vecin(nz,ny,nx), tmp(nz,ny,nx)
      common/io/inp, iout
      call ebc(tmp,hz,vecin,nz,nz,ny*nx)
      do 70 j=1,nx
         do 80 k=1,ny
            do 90 l=1,nz
               vecout(l,k,j) = vecout(l,k,j) + tmp(l,k,j) 
 90         continue
 80      continue   
 70   continue
      do 200 i=1,nx
         call ebc(tmp(1,1,i),vecin(1,1,i),hy,nz,ny,ny)
 200  continue   
      do 400 j=1,nx
         do 500 k=1,ny
            do 600 l=1,nz
               vecout(l,k,j) = vecout(l,k,j) + tmp(l,k,j) 
 600        continue
 500     continue   
 400  continue   
      call ebc(tmp,vecin,hx,nz*ny,nx,nx)
      do 900 j=1,nx
         do 1000 k=1,ny
            do 1100 l=1,nz
               vecout(l,k,j) = vecout(l,k,j) + tmp(l,k,j) 
 1100       continue
 1000    continue   
 900  continue   
      return
      end       
