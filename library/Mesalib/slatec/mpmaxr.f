*deck mpmaxr
      subroutine mpmaxr (x)
c***begin prologue  mpmaxr
c***subsidiary
c***purpose  subsidiary to dqdota and dqdoti
c***library   slatec
c***type      all (mpmaxr-a)
c***author  (unknown)
c***description
c
c  sets x to the largest possible positive 'mp' number.
c
c  the argument x(*) is an integer arrays of size 30.  see the comments
c  in the routine mpblas for the reason for this choice.
c
c***see also  dqdota, dqdoti, mpblas
c***routines called  mpchk
c***common blocks    mpcom
c***revision history  (yymmdd)
c   791001  date written
c   ??????  modified for use with blas.  blank common changed to named
c           common.  r given dimension 12.
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   930124  increased array size in mpcon for sun -r8.  (rwc)
c***end prologue  mpmaxr
      common /mpcom/ b, t, m, lun, mxr, r(30)
      integer b, t, r, x(*)
c***first executable statement  mpmaxr
      call mpchk (1, 4)
      it = b - 1
c set fraction digits to b-1
      do 10 i = 1, t
   10 x(i+2) = it
c set sign and exponent
      x(1) = 1
      x(2) = m
      return
      end
