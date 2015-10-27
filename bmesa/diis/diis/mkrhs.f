*deck mkrhs.f
c***begin prologue     mkrhs
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
c***end prologue       mkrhs
      subroutine mkrhs(rhs,energy,p,x,wt,v,type,n,npt,nrhs)
      implicit integer (a-z)
      real*8 rhs
      real*8 p, x, wt, v, energy, k, tmp, sdot
      character*(*) type
      dimension p(npt,0:n-1), x(npt), wt(npt), v(npt), rhs(n,nrhs) 
      common/io/inp, iout
      k=sqrt(2.d0*energy)
      if(type.eq.'exponential') then
         do 10 i=1,npt
            v(i) =  - wt(i)*exp(-x(i))*sin(k*x(i))
 10      continue
      elseif(type.eq.'well') then
         do 20 i=1,npt
            v(i) =  - wt(i)*1.0d0*sin(k*x(i))
 20      continue
      endif
      call rzero(rhs,n*nrhs)
      do 30 i=1,n
         tmp=sdot(npt,p(1,i-1),1,v,1)
         rhs(i,1)=tmp
 30   continue   
      return
      end       
