*deck dcbrt
      double precision function dcbrt (x)
c***begin prologue  dcbrt
c***purpose  compute the cube root.
c***library   slatec (fnlib)
c***category  c2
c***type      double precision (cbrt-s, dcbrt-d, ccbrt-c)
c***keywords  cube root, elementary functions, fnlib, roots
c***author  fullerton, w., (lanl)
c***description
c
c dcbrt(x) calculates the double precision cube root for
c double precision argument x.
c
c***references  (none)
c***routines called  d1mach, d9pak, d9upak
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  dcbrt
      double precision x, cbrt2(5), y, cbrtsq,  d9pak, d1mach
      save cbrt2, niter
      data cbrt2(1) / 0.6299605249 4743658238 3605303639 11 d0 /
      data cbrt2(2) / 0.7937005259 8409973737 5852819636 15 d0 /
      data cbrt2(3) / 1.0 d0 /
      data cbrt2(4) / 1.2599210498 9487316476 7210607278 23 d0 /
      data cbrt2(5) / 1.5874010519 6819947475 1705639272 31 d0 /
      data niter / 0 /
c***first executable statement  dcbrt
      if (niter.eq.0) niter = 1.443*log(-.106*log(0.1*real(d1mach(3)))
     1  ) + 1.0
c
      dcbrt = 0.d0
      if (x.eq.0.d0) return
c
      call d9upak (abs(x), y, n)
      ixpnt = n/3
      irem = n - 3*ixpnt + 3
c
c the approximation below is a generalized chebyshev series converted
c to polynomial form.  the approx is nearly best in the sense of
c relative error with 4.085 digits accuracy.
c
      z = y
      dcbrt = .439581e0 + z*(.928549e0 + z*(-.512653e0 + z*.144586e0))
c
      do 10 iter=1,niter
        cbrtsq = dcbrt*dcbrt
        dcbrt = dcbrt + (y-dcbrt*cbrtsq)/(3.d0*cbrtsq)
 10   continue
c
      dcbrt = d9pak (cbrt2(irem)*sign(dcbrt,x), ixpnt)
      return
c
      end
