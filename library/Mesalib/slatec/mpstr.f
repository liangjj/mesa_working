*deck mpstr
      subroutine mpstr (x, y)
c***begin prologue  mpstr
c***subsidiary
c***purpose  subsidiary to dqdota and dqdoti
c***library   slatec
c***type      all (mpstr-a)
c***author  (unknown)
c***description
c
c  sets y = x for 'mp' x and y.
c
c  the arguments x(*) and y(*) are integer arrays of size 30.  see the
c  comments in the routine mpblas for the reason for this choice.
c
c***see also  dqdota, dqdoti, mpblas
c***routines called  (none)
c***common blocks    mpcom
c***revision history  (yymmdd)
c   791001  date written
c   ??????  modified for use with blas.  blank common changed to named
c           common.  r given dimension 12.
c   890206  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   930124  increased array size in mpcon for sun -r8.  (rwc)
c***end prologue  mpstr
      common /mpcom/ b, t, m, lun, mxr, r(30)
      integer b, t, r, x(*), y(*)
c***first executable statement  mpstr
      do 10 i = 1, t+2
         y(i) = x(i)
   10 continue
      return
      end
