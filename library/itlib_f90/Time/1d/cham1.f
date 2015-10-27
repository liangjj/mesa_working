*deck cham1.f
c***begin prologue     cham1
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
c***                   the index array is computed for the full hamiltonian
c***                   which is why its column length is m not n.
c***                   
c***references         
c
c***routines called    
c***end prologue       cham1
      subroutine cham1(h,hx,ht,vxt,ind,nx,nt,n,m,prn)
      implicit integer (a-z)
      complex*16 h, eye
      real*8 hx, ht, vxt
      character*80 title 
      logical prn
      dimension h(n,n), hx(nx,nx), ht(nt,nt), vxt(n), ind(m,2)
      common/io/inp, iout
      data eye/(0.d0,1.d0)/
      call czero(h,n*n)
      do 10 i=1,n
         ti=ind(i,1) 
         xi=ind(i,2)
         do 20 j=1,n
            tj=ind(j,1)
            xj=ind(j,2)
            dltitj=0
            dlxixj=0
            if(ti.eq.tj) then
               dltitj=1
            endif 
            if(xi.eq.xj) then
               dlxixj=1
            endif 
            h(i,j) = h(i,j) - hx(xi,xj)*dltitj + eye*ht(ti,tj)*dlxixj
 20      continue
 10   continue   
      do 30 i=1,n
         h(i,i) = h(i,i) - vxt(i)
 30   continue   
      if(prn) then 
         title='complex matrix'
         call prntcm(title,h,n,n,n,n,iout) 
      endif
      return
      end       

