*deck tstsol.f
c***begin prologue     tstsol
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            test time-dependent wavefunction by expansion.
c***                   
c***references         
c
c***routines called    
c***end prologue       tstsol
      subroutine tstsol(p,q,rhs,n,npts,type)
      implicit integer (a-z)
      real*8 p, q
      complex*16 rhs, exact, eye, temp
      character*80 title
      character*(*) type
      dimension p(npts,0:n-1), q(npts), rhs(n) 
      common/io/inp, iout
      data eye/(0.d0,1.d0)/
      if(type.eq.'1') then
         do 10 i=1,npts
            exact = q(i) + eye*( exp(eye*q(i)) - 1.d0 )
            temp=(0.d0,0.d0)
            do 20 j=1,n
               temp = temp + p(i,j-1)*rhs(j)
 20         continue
            write(iout,1) q(i), exact, temp
 10      continue   
      elseif(type.eq.'t') then
         do 30 i=1,npts
            exact = exp(-.5d0*eye*q(i)*q(i)) - 1.d0
            temp=(0.d0,0.d0)
            do 40 j=1,n
               temp = temp + p(i,j-1)*rhs(j)
 40         continue
            write(iout,1) q(i), exact, temp
 30      continue   
      endif
      return
 1    format(1x,'time = ',e15.8,/,1x,
     1          'exact value      = ',e15.8,1x,e15.8,/,1x,
     2          'calculated value = ',e15.8,1x,e15.8)
      end       
