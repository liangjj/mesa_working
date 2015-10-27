*deck hamfl2.f
c***begin prologue     hamfl2
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
c***end prologue       hamfl2
      subroutine hamfl2(h,hx,hy,vxy,ind,nx,ny,n,prn)
      implicit integer (a-z)
      real*8 h, hx, hy, vxy
      logical prn
      dimension h(n,n), ind(ny,nx), hx(nx,nx), hy(ny,ny), vxy(n)
      common/io/inp, iout
      do 10 xi=1,nx
         do 20 yi=1,ny
            yixi=ind(yi,xi)
            do 30 yj=1,ny
               yjxi=ind(yj,xi)
               h(yixi,yjxi) = h(yixi,yjxi) + hy(yi,yj)
 30         continue
            do 40 xj=1,nx
               yixj=ind(yi,xj)   
               h(yixi,yixj) = h(yixi,yixj) + hx(xi,xj)
 40            continue   
 20      continue
 10   continue   
      do 100 xi=1,nx
         do 200 yi=1,ny
            yixi=ind(yi,xi)
            h(yixi,yixi) = h(yixi,yixi) + vxy(yixi)   
 200     continue
 100  continue   
      return
      end       


