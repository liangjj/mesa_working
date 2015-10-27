*deck rham1.f
c***begin prologue     rham1
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form all or a sub-block of the two dimensional 
c***                   real hamiltonian.
c***                   
c***                   in setting up this sub-block of the hamiltonian
c***                   there is the explicit condition that it not cross
c***                   channel boundaries.
c***                   the index array is computed for the full hamiltonian
c***                   which is why its column length is m not n.
c***                   
c***references         
c
c***routines called    
c***end prologue       rham1
      subroutine rham1(h,hx,vx,ind,nx,n,m,prn)
      implicit integer (a-z)
      real*8 h, hx, vx
      character*80 title 
      logical prn
      dimension h(n,n), hx(nx,nx), vx(n), ind(m,2)
      common/io/inp, iout
      call rzero(h,n*n)
      do 10 i=1,n
         xi=ind(i,1)
         do 20 j=1,n
            xj=ind(j,1)
            h(i,j) = h(i,j) + hx(xi,xj) 
 20      continue
 10   continue   
      do 30 i=1,n
         h(i,i) = h(i,i) + vx(i)
 30   continue   
      if(prn) then 
         title='matrix'
         call prntrm(title,h,n,n,n,n,iout) 
      endif
      return
      end       

