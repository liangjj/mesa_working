*deck vhmohy.f
c***begin prologue     vhmohy
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           potential, hyperspherical
c***author             schneider, barry (nsf)
c***source             
c***purpose            harmonic oscillator potential in hyperspherical coordinates
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       vhmohy
      subroutine vhmohy(v,phi,r,mass,omega,n)
      implicit integer (a-z)
      real*8 v, r, phi
      real*8 r1, r2, r12, omega, mass
      real*8 fpkey, ci, si
      dimension omega(2,2), mass(2,2)
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
            v(i,j) = v(i,j) + .5d0*mass(1,1)*omega(1,1)*omega(1,1)*r1*r1
     1                      + .5d0*mass(2,2)*omega(2,2)*omega(2,2)*r2*r2
     2                      + .5d0*mass(1,2)*omega(1,2)*omega(1,2)
     3                                                 *r12*r12
 20      continue
 10   continue
      return
      end       














