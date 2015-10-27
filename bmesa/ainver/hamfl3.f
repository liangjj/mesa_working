*deck hamfl3.f
c***begin prologue     hamfl3
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form three dimensional hamiltonian.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       hamfl3
      subroutine hamfl3(h,hx,hy,hz,vxyz,ind,nx,ny,nz,n)
      implicit integer (a-z)
      real*8 h, hx, hy, hz, vxyz
      dimension h(n,n), ind(nz,ny,nx), hx(nx,nx), hy(ny,ny), hz(nz,nz)
      dimension vxyz(n)
      common/io/inp, iout
      do 10 xi=1,nx
         do 20 yi=1,ny
            do 30 zi=1,nz
               ziyixi=ind(zi,yi,xi)
               do 40 zj=1,nz
                  zjyixi=ind(zj,yi,xi)
                  h(ziyixi,zjyixi) = h(ziyixi,zjyixi) + hz(zi,zj)
 40            continue
               do 50 yj=1,ny
                  ziyjxi=ind(zi,yj,xi)
                  h(ziyixi,ziyjxi) = h(ziyixi,ziyjxi) + hy(yi,yj)
 50            continue
               do 60 xj=1,nx
                  ziyixj=ind(zi,yi,xj)
                  h(ziyixi,ziyixj) = h(ziyixi,ziyixj) + hx(xi,xj)
 60            continue   
 30         continue   
 20      continue
 10   continue   
      do 100 xi=1,nx
         do 200 yi=1,ny
            do 300 zi=1,nz
               ziyixi=ind(zi,yi,xi)
               h(ziyixi,ziyixi) = h(ziyixi,ziyixi) + vxyz(ziyixi)   
 300        continue   
 200     continue
 100  continue   
      return
      end       


