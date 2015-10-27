*deck mpcmd
      subroutine mpcmd (x, dz)
c***begin prologue  mpcmd
c***subsidiary
c***purpose  subsidiary to dqdota and dqdoti
c***library   slatec
c***type      all (mpcmd-a)
c***author  (unknown)
c***description
c
c  converts multiple-precision x to double-precision dz. assumes
c  x is in allowable range for double-precision numbers. there is
c  some loss of accuracy if the exponent is large.
c
c  the argument x(*) is integer array of size 30.  see the comments in
c  the routine mpblas for the reason for this choice.
c
c***see also  dqdota, dqdoti, mpblas
c***routines called  mpchk, mperr
c***common blocks    mpcom
c***revision history  (yymmdd)
c   791001  date written
c   ??????  modified for use with blas.  blank common changed to named
c           common.  r given dimension 12.
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   930124  increased array size in mpcon for sun -r8.  (rwc)
c***end prologue  mpcmd
      double precision db, dz, dz2
      common /mpcom/ b, t, m, lun, mxr, r(30)
      integer b, t, r, x(*), tm
c***first executable statement  mpcmd
      call mpchk (1, 4)
      dz = 0d0
      if (x(1).eq.0) return
      db = dble(b)
      do 10 i = 1, t
      dz = db*dz + dble(x(i+2))
      tm = i
c check if full double-precision accuracy attained
      dz2 = dz + 1d0
c test below not always equivalent to - if (dz2.le.dz) go to 20,
c for example on cyber 76.
      if ((dz2-dz).le.0d0) go to 20
   10 continue
c now allow for exponent
   20 dz = dz*(db**(x(2)-tm))
c check reasonableness of result.
      if (dz.le.0d0) go to 30
c lhs should be .le. 0.5 but allow for some error in log
      if (abs(dble(x(2))-(log(dz)/
     1    log(dble(b))+0.5d0)).gt.0.6d0) go to 30
      if (x(1).lt.0) dz = -dz
      return
c following message indicates that x is too large or small -
c try using mpcmde instead.
   30 write (lun, 40)
   40 format (' *** floating-point over/under-flow in mpcmd ***')
      call mperr
      return
      end
