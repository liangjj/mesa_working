*deck rham3.f
c***begin prologue     rham3
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form all or a sub-block of the three dimensional 
c***                   
c***                   the index array is computed for the full hamiltonian
c***                   which is why its column length is m not n.
c***references         
c
c***routines called    
c***end prologue       rham3
      subroutine rham3(h,hx,hy,hz,vxyz,ind,nx,ny,nz,n,m,prn)
      implicit integer (a-z)
      real*8 h, hx, hy, hz, vxyz
      character*80 title 
      logical prn
      dimension h(n,n), hx(nx,nx), hy(ny,ny), hz(nz,nz)
      dimension vxyz(n), ind(m,3)
      common/io/inp, iout
      call rzero(h,n*n)
      do 10 i=1,n
         zi=ind(i,1)
         yi=ind(i,2) 
	 xi=ind(i,3)
         do 20 j=1,n
            zj=ind(j,1)
            yj=ind(j,2)
            xj=ind(j,3) 
            dlxixj=0
            dlyiyj=0
            dlzizj=0
            if(zi.eq.zj) then
               dlzizj=1
            endif 
            if(xi.eq.xj) then
               dlxixj=1
            endif 
            if(yi.eq.yj) then
               dlyiyj=1
            endif 
            h(i,j) = h(i,j) + hx(xi,xj)*dlzizj*dlyiyj 
     1                      + hy(yi,yj)*dlzizj*dlxixj 
     2                      + hz(zi,zj)*dlxixj*dlyiyj
 20      continue
 10   continue   
      do 30 i=1,n
         h(i,i) = h(i,i) + vxyz(i)
 30   continue   
      if(prn) then 
         title='real matrix'
         call prntcm(title,h,n,n,n,n,iout) 
      endif
      return
      end       

