*deck c9ln2r
      complex function c9ln2r (z)
c***begin prologue  c9ln2r
c***subsidiary
c***purpose  evaluate log(1+z) from second order relative accuracy so
c            that  log(1+z) = z - z**2/2 + z**3*c9ln2r(z).
c***library   slatec (fnlib)
c***category  c4b
c***type      complex (r9ln2r-s, d9ln2r-d, c9ln2r-c)
c***keywords  elementary functions, fnlib, logarithm, second order
c***author  fullerton, w., (lanl)
c***description
c
c evaluate  log(1+z)  from 2-nd order with relative error accuracy so
c that     log(1+z) = z - z**2/2 + z**3*c9ln2r(z).
c
c now  log(1+z) = 0.5*log(1+2*x+abs(z)**2) + i*carg(1+z),
c where x = real(z)  and  y = aimag(z).
c we find
c     z**3 * c9ln2r(z) = -x*abs(z)**2 - 0.25*abs(z)**4
c        + (2*x+abs(z)**2)**3 * r9ln2r(2*x+abs(z)**2)
c        + i * (carg(1+z) + (x-1)*y)
c the imaginary part must be evaluated carefully as
c     (atan(y/(1+x)) - y/(1+x)) + y/(1+x) - (1-x)*y
c       = (y/(1+x))**3 * r9atn1(y/(1+x)) + x**2*y/(1+x)
c
c now we divide through by z**3 carefully.  write
c     1/z**3 = (x-i*y)/abs(z)**3 * (1/abs(z)**3)
c then   c9ln2r(z) = ((x-i*y)/abs(z))**3 * (-x/abs(z) - abs(z)/4
c        + 0.5*((2*x+abs(z)**2)/abs(z))**3 * r9ln2r(2*x+abs(z)**2)
c        + i*y/(abs(z)*(1+x)) * ((x/abs(z))**2 +
c          + (y/(abs(z)*(1+x)))**2 * r9atn1(y/(1+x)) ) )
c
c if we let  xz = x/abs(z)  and  yz = y/abs(z)  we may write
c     c9ln2r(z) = (xz-i*yz)**3 * (-xz - abs(z)/4
c        + 0.5*(2*xz+abs(z))**3 * r9ln2r(2*x+abs(z)**2)
c        + i*yz/(1+x) * (xz**2 + (yz/(1+x))**2*r9atn1(y/(1+x)) ))
c
c***references  (none)
c***routines called  r9atn1, r9ln2r
c***revision history  (yymmdd)
c   780401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900720  routine changed from user-callable to subsidiary.  (wrb)
c***end prologue  c9ln2r
      complex z
c***first executable statement  c9ln2r
      x = real (z)
      y = aimag (z)
c
      cabsz = abs(z)
      if (cabsz.gt.0.8125) go to 20
c
      c9ln2r = cmplx (1.0/3.0, 0.0)
      if (cabsz.eq.0.0) return
c
      xz = x/cabsz
      yz = y/cabsz
c
      arg = 2.0*xz + cabsz
      rpart = 0.5*arg**3*r9ln2r(cabsz*arg) - xz - 0.25*cabsz
      y1x = yz/(1.0+x)
      aipart = y1x * (xz**2 + y1x**2*r9atn1(cabsz*y1x) )
c
      c9ln2r = cmplx(xz,-yz)**3 * cmplx(rpart,aipart)
      return
c
 20   c9ln2r = (log(1.0+z) - z*(1.0-0.5*z)) / z**3
      return
c
      end
