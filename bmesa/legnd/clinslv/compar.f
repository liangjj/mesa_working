*deck compar.f
c***begin prologue     compar
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            compare numerical and exact solutions
c***                   for a well interaction with either zero function
c***                   or zero derivative boundary conditions at right end.            
c***                   
c***references         
c
c***routines called    
c***end prologue       compar
      subroutine compar(q,wt,p,rhs,energy,rtbc,n)
      implicit integer (a-z)
      real*8 q, wt, p, energy, k, fac1, fac2
      complex*16 rhs, exact, approx, eye, cdotu
      dimension q(n), wt(n), p(n,n), rhs(n) 
      common/io/inp, iout
      data eye/(0.d0,1.d0)/
      k=sqrt(2.d0*energy)
      if(rtbc.eq.0) then
         write(iout,*) ' zero function boundary condition'
         fac1=sin(k)
         fac2=2.d0/(k*k*fac1)      
         do 10 i=1,n
            exact = fac2*(sin(k*(q(i)-1.d0)) - sin(k*q(i)) + fac1 )
            approx = rhs(i)*p(i,i)   
            write(iout,1) q(i), exact, approx
   10    continue
      else
         write(iout,*) ' zero derivative boundary condition'
         fac1=cos(k)
         fac2=1.d0/(k*k*fac1)
         fac1=1.d0/(k*k)      
         do 30 i=1,n
            exact = 2.d0*(fac1 - fac2*cos(k*(q(i)-1.d0)))
            approx = rhs(i)*p(i,i)
            write(iout,1) q(i), exact, approx
 30      continue
      endif    
 1    format(1x,'x = ',e15.8,/,5x,'exact      = (',e15.8,',',e15.8,')',
     1       /,5x,'approximate = (',e15.8,',',e15.8,')')
      return
      end       
