*deck dgamr
      double precision function dgamr (x)
c***begin prologue  dgamr
c***purpose  compute the reciprocal of the gamma function.
c***library   slatec (fnlib)
c***category  c7a
c***type      double precision (gamr-s, dgamr-d, cgamr-c)
c***keywords  fnlib, reciprocal gamma function, special functions
c***author  fullerton, w., (lanl)
c***description
c
c dgamr(x) calculates the double precision reciprocal of the
c complete gamma function for double precision argument x.
c
c***references  (none)
c***routines called  dgamma, dlgams, xerclr, xgetf, xsetf
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900727  added external statement.  (wrb)
c***end prologue  dgamr
      double precision x, alngx, sgngx, dgamma
      external dgamma
c***first executable statement  dgamr
      dgamr = 0.0d0
      if (x.le.0.0d0 .and. aint(x).eq.x) return
c
      call xgetf (irold)
      call xsetf (1)
      if (abs(x).gt.10.0d0) go to 10
      dgamr = 1.0d0/dgamma(x)
      call xerclr
      call xsetf (irold)
      return
c
 10   call dlgams (x, alngx, sgngx)
      call xerclr
      call xsetf (irold)
      dgamr = sgngx * exp(-alngx)
      return
c
      end
