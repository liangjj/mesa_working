*deck mkpsi0.f
c***begin prologue     mkpsi0
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***references         
c
c***routines called    
c***end prologue       mkpsi0
      subroutine mkpsi0(rhs,energy,p,x,wt,v,type,n,npt,nrhs)
      implicit integer (a-z)
      complex*16 rhs, v, eye
      real*8 p, x, wt, energy, k
      character*(*) type
      dimension p(npt,0:n-1), x(npt), wt(npt), v(npt), rhs(n,nrhs) 
      common/io/inp, iout
      data eye /(0.d0,1.d0)/
      k=sqrt(2.d0*energy)
      if(type.eq.'exponential') then
         do 10 i=1,npt
            v(i) =  - wt(i)*exp(-x(i))*sin(k*x(i))
 10      continue
      elseif(type.eq.'well') then
         do 20 i=1,npt
            v(i) =  - wt(i)*1.0d0*sin(k*x(i))
 20      continue
       elseif(type.eq.'complex-well') then
         do 30 i=1,n
            v(i) = ( -1.d0 - eye*.2d0 )*wt(i)*sin(k*x(i))
 30      continue            
      endif
      call czero(rhs,n*nrhs)
      do 40 i=1,n
         do 50 j=1,npt
            rhs(i,1) = rhs(i,1) + p(j,i-1)*v(j) 
 50      continue   
 40   continue   
      return
      end       
