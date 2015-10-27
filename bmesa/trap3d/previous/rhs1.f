*deck rhs1.f
c***begin prologue     rhs1
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            right hand side of inhomogeneous time-dependent
c***                   hamiltonian.
c***                   
c***description        simple time-dependent hamiltonian using polynomial
c***                   basis.  
c***references         
c
c***routines called    
c***end prologue       rhs1
      subroutine rhs1(p,q,wt,rhs,n,npts)
      implicit integer (a-z)
      real*8 p, q, wt
      complex*16  rhs
      dimension q(npts), p(npts,0:n-1), wt(npts), rhs(n) 
      common/io/inp, iout
      call czero(rhs,n)
      do 10 i=1,n
         do 20 k=1,npts
            rhs(i) = rhs(i) +wt(k)*q(k)*p(k,i-1)
 20      continue
 10   continue   
      return
      end       
