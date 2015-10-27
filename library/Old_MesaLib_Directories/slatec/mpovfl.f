*deck mpovfl
      subroutine mpovfl (x)
c***begin prologue  mpovfl
c***subsidiary
c***purpose  subsidiary to dqdota and dqdoti
c***library   slatec
c***type      all (mpovfl-a)
c***author  (unknown)
c***description
c
c  called on multiple-precision overflow, i.e. when the
c  exponent of 'mp' number x would exceed m.  at present execution is
c  terminated with an error message after calling mpmaxr(x), but it
c  would be possible to return, possibly updating a counter and
c  terminating execution after a preset number of overflows.  action
c  could easily be determined by a flag in labelled common.
c
c  the argument x(*) is an integer array of size 30.  see the comments
c  in the routine mpblas for the reason for this choice.
c
c***see also  dqdota, dqdoti, mpblas
c***routines called  mpchk, mperr, mpmaxr
c***common blocks    mpcom
c***revision history  (yymmdd)
c   791001  date written
c   ??????  modified for use with blas.  blank common changed to named
c           common.  r given dimension 12.
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   930124  increased array size in mpcon for sun -r8.  (rwc)
c***end prologue  mpovfl
      common /mpcom/ b, t, m, lun, mxr, r(30)
      integer b, t, r, x(*)
c***first executable statement  mpovfl
      call mpchk (1, 4)
c set x to largest possible positive number
      call mpmaxr (x)
      write (lun, 10)
   10 format (' *** call to mpovfl, mp overflow occurred ***')
c terminate execution by calling mperr
      call mperr
      return
      end
