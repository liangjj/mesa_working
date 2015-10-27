*deck cham3.f
c***begin prologue     cham3
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form all or a sub-block of the two dimensional 
c***                   complex hamiltonian.
c***                   
c***                   in setting up this sub-block of the hamiltonian
c***                   there is the explicit condition that it not cross
c***                   channel boundaries.
c***
c***                   the index array is computed for the full hamiltonian
c***                   which is why its column length is m not n.
c***references         
c
c***routines called    
c***end prologue       cham3
      subroutine cham3(h,hx,hy,hz,ht,vxyzt,ind,nx,ny,nz,nt,n,m,prn)
      implicit integer (a-z)
      complex*16 h, eye
      real*8 hx, hy, hz, ht, vxyzt
      character*80 title 
      logical prn
      dimension h(n,n), hx(nx,nx), hy(ny,ny), hz(nz,nz)
      dimension ht(nt,nt), vxyzt(n), ind(m,4)
      dimension hdum(2*n,2*n)
      common/io/inp, iout
      data eye/(0.d0,1.d0)/
      call czero(h,n*n)
      do 10 i=1,n
         xi=ind(i,4) 
         yi=ind(i,3) 
         zi=ind(i,2) 
         ti=ind(i,1)
         do 20 j=1,n
            xj=ind(j,4) 
            yj=ind(j,3) 
            zj=ind(j,2) 
            tj=ind(j,1)
            dlxixj=0
            dlyiyj=0
            dlzizj=0
            dltitj=0
            if(ti.eq.tj) then
               dltitj=1
            endif 
            if(xi.eq.xj) then
               dlxixj=1
            endif 
            if(yi.eq.yj) then
               dlyiyj=1
            endif 
            if(zi.eq.zj) then
               dlzizj=1
            endif 
            h(i,j) = h(i,j) - hx(xi,xj)*dltitj*dlyiyj*dlzizj 
     1                      - hy(yi,yj)*dltitj*dlxixj*dlzizj
     2                      - hz(zi,zj)*dltitj*dlxixj*dlyiyj 
     2                      + eye*ht(ti,tj)*dlxixj*dlyiyj*dlzizj
 20      continue
 10   continue   
      do 30 i=1,n
         h(i,i) = h(i,i) - vxyzt(i)
 30   continue   
      if(prn) then 
         title='complex matrix'
         call prntcm(title,h,n,n,n,n,iout) 
      endif
      return
      end       

