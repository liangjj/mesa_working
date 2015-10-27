*deck c0lgmc
      complex function c0lgmc (z)
c***begin prologue  c0lgmc
c***purpose  evaluate (z+0.5)*log((z+1.)/z) - 1.0 with relative
c            accuracy.
c***library   slatec (fnlib)
c***category  c7a
c***type      complex (c0lgmc-c)
c***keywords  fnlib, gamma function, special functions
c***author  fullerton, w., (lanl)
c***description
c
c evaluate  (z+0.5)*log((z+1.0)/z) - 1.0  with relative error accuracy
c let q = 1.0/z so that
c     (z+0.5)*log(1+1/z) - 1 = (z+0.5)*(log(1+q) - q + q*q/2) - q*q/4
c        = (z+0.5)*q**3*c9ln2r(q) - q**2/4,
c where  c9ln2r  is (log(1+q) - q + 0.5*q**2) / q**3.
c
c***references  (none)
c***routines called  c9ln2r, r1mach
c***revision history  (yymmdd)
c   780401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  c0lgmc
      complex z, q, c9ln2r
      save rbig
      data rbig / 0.0 /
c***first executable statement  c0lgmc
      if (rbig.eq.0.0) rbig = 1.0/r1mach(3)
c
      cabsz = abs(z)
      if (cabsz.gt.rbig) c0lgmc = -(z+0.5)*log(z) - z
      if (cabsz.gt.rbig) return
c
      q = 1.0/z
      if (cabsz.le.1.23) c0lgmc = (z+0.5)*log(1.0+q) - 1.0
      if (cabsz.gt.1.23) c0lgmc = ((1.+.5*q)*c9ln2r(q) - .25) * q**2
c
      return
      end
