*deck rham2.f
c***begin prologue     rham2
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form all or a sub-block of the two dimensional 
c***                   hamiltonian.
c***                   
c***                   the index array is computed for the full hamiltonian
c***                   which is why its column length is m not n.
c***references         
c
c***routines called    
c***end prologue       rham2
      subroutine rham2(h,hx,hy,vxy,ind,nx,ny,n,m,prn)
      implicit integer (a-z)
      real*8 h, hx, hy, vxy
      character*80 title 
      logical prn
      dimension h(n,n), hx(nx,nx), hy(ny,ny)
      dimension vxy(n), ind(m,2)
      common/io/inp, iout
      call rzero(h,n*n)
      do 10 i=1,n
         yi=ind(i,1) 
	 xi=ind(i,2)
         do 20 j=1,n
            yj=ind(j,1)
            xj=ind(j,2) 
            dlxixj=0
            dlyiyj=0
            if(xi.eq.xj) then
               dlxixj=1
            endif 
            if(yi.eq.yj) then
               dlyiyj=1
            endif 
            h(i,j) = h(i,j) + hx(xi,xj)*dlyiyj 
     1                      + hy(yi,yj)*dlxixj
 20      continue
 10   continue   
      do 30 i=1,n
         h(i,i) = h(i,i) + vxy(i)
 30   continue   
      if(prn) then 
         title='real matrix'
         call prntrm(title,h,n,n,n,n,iout) 
      endif
      return
      end       

