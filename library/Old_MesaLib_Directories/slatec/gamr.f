*deck gamr
      function gamr (x)
c***begin prologue  gamr
c***purpose  compute the reciprocal of the gamma function.
c***library   slatec (fnlib)
c***category  c7a
c***type      single precision (gamr-s, dgamr-d, cgamr-c)
c***keywords  fnlib, reciprocal gamma function, special functions
c***author  fullerton, w., (lanl)
c***description
c
c gamr is a single precision function that evaluates the reciprocal
c of the gamma function for single precision argument x.
c
c***references  (none)
c***routines called  algams, gamma, xerclr, xgetf, xsetf
c***revision history  (yymmdd)
c   770701  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900727  added external statement.  (wrb)
c***end prologue  gamr
      external gamma
c***first executable statement  gamr
      gamr = 0.0
      if (x.le.0.0 .and. aint(x).eq.x) return
c
      call xgetf (irold)
      call xsetf (1)
      if (abs(x).gt.10.0) go to 10
      gamr = 1.0/gamma(x)
      call xerclr
      call xsetf (irold)
      return
c
 10   call algams (x, alngx, sgngx)
      call xerclr
      call xsetf (irold)
      gamr = sgngx * exp(-alngx)
      return
c
      end
