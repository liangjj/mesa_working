*deck vexphy.f
c***begin prologue     vexphy
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           potential, hyperspherical
c***author             schneider, barry (nsf)
c***source             
c***purpose            exponential potential in hyperspherical coordinates
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       vexphy
      subroutine vexphy(v,phi,r,d,a,n)
      implicit integer (a-z)
      real*8 v, r, phi
      real*8 r1, r2, r12, d, a, ci, si
      dimension d(2,2), a(2,2)
      dimension n(2)
      dimension v(n(1),n(2)), r(n(2)), phi(n(1))
      common/io/inp, iout
      do 10 i=1,n(1)
         ci=cos(phi(i))
         si=sin(phi(i))
         do 20 j=1,n(2)
            r1=r(j)*ci
            r2=r(j)*si
            r12=abs(r2-r1)
            v(i,j) = v(i,j) + d(1,1)*exp(-a(1,1)*r1)
     1                      + d(2,2)*exp(-a(2,2)*r1)
     2                      + d(1,2)*exp(-a(1,2)*r12)
 20         continue
 10      continue    
      return
      end       














