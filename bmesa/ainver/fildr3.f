*deck fildr3.f
c***begin prologue     fildr3
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***routines called    
c***end prologue       fildr3
      subroutine fildr3(driver,vxyz,psi0,nx,ny,nz)
      implicit integer (a-z)
      real*8 vxyz, driver, psi0
      dimension driver(nz,ny,nx), psi0(nz,ny,nx)
      dimension vxyz(nz,ny,nx)
      common/io/inp, iout
      do 20 j=1,nx
         do 30 k=1,ny
            do 40 l=1,nz
               driver(l,k,j) = - vxyz(l,k,j)*psi0(l,k,j)
 40         continue
 30      continue
 20   continue   
      return
      end       
