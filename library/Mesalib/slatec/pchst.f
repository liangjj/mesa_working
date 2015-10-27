*deck pchst
      real function pchst (arg1, arg2)
c***begin prologue  pchst
c***subsidiary
c***purpose  pchip sign-testing routine
c***library   slatec (pchip)
c***type      single precision (pchst-s, dpchst-d)
c***author  fritsch, f. n., (llnl)
c***description
c
c         pchst:  pchip sign-testing routine.
c
c     returns:
c        -1. if arg1 and arg2 are of opposite sign.
c         0. if either argument is zero.
c        +1. if arg1 and arg2 are of the same sign.
c
c     the object is to do this without multiplying arg1*arg2, to avoid
c     possible over/underflow problems.
c
c  fortran intrinsics used:  sign.
c
c***see also  pchce, pchci, pchcs, pchim
c***routines called  (none)
c***revision history  (yymmdd)
c   811103  date written
c   820805  converted to slatec library version.
c   870813  minor cosmetic changes.
c   890411  added save statements (vers. 3.2).
c   890411  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910408  updated author and date written sections in prologue.  (wrb)
c   930503  improved purpose.  (fnf)
c***end prologue  pchst
c
c**end
c
c  declare arguments.
c
      real  arg1, arg2
c
c  declare local variables.
c
      real  one, zero
      save zero, one
      data  zero /0./,  one /1./
c
c  perform the test.
c
c***first executable statement  pchst
      pchst = sign(one,arg1) * sign(one,arg2)
      if ((arg1.eq.zero) .or. (arg2.eq.zero))  pchst = zero
c
      return
c------------- last line of pchst follows ------------------------------
      end
