*deck cbrt
      function cbrt (x)
c***begin prologue  cbrt
c***purpose  compute the cube root.
c***library   slatec (fnlib)
c***category  c2
c***type      single precision (cbrt-s, dcbrt-d, ccbrt-c)
c***keywords  cube root, elementary functions, fnlib, roots
c***author  fullerton, w., (lanl)
c***description
c
c cbrt(x) calculates the cube root of x.
c
c***references  (none)
c***routines called  r1mach, r9pak, r9upak
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  cbrt
      dimension cbrt2(5)
      save cbrt2, niter
      data cbrt2(1) / 0.6299605249 4743658e0 /
      data cbrt2(2) / 0.7937005259 8409974e0 /
      data cbrt2(3) / 1.0e0 /
      data cbrt2(4) / 1.2599210498 9487316e0 /
      data cbrt2(5) / 1.5874010519 6819947e0 /
      data niter / 0 /
c***first executable statement  cbrt
      if (niter.eq.0) niter = 1.443*log(-.106*log(0.1*r1mach(3))) + 1.
c
      cbrt = 0.0
      if (x.eq.0.) return
c
      call r9upak (abs(x), y, n)
      ixpnt = n/3
      irem = n - 3*ixpnt + 3
c
c the approximation below is a generalized chebyshev series converted
c to polynomial form.  the approx is nearly best in the sense of
c relative error with 4.085 digits accuracy.
c
      cbrt = .439581e0 + y*(.928549e0 + y*(-.512653e0 + y*.144586e0))
c
      do 10 iter=1,niter
        cbrtsq = cbrt*cbrt
        cbrt = cbrt + (y-cbrt*cbrtsq)/(3.0*cbrtsq)
 10   continue
c
      cbrt = r9pak (cbrt2(irem)*sign(cbrt,x), ixpnt)
      return
c
      end
