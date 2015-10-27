*deck mpunfl
      subroutine mpunfl (x)
c***begin prologue  mpunfl
c***subsidiary
c***purpose  subsidiary to dqdota and dqdoti
c***library   slatec
c***type      all (mpunfl-a)
c***author  (unknown)
c***description
c
c called on multiple-precision underflow, i.e.  when the
c exponent of 'mp' number x would be less than -m.
c
c***see also  dqdota, dqdoti
c***routines called  mpchk
c***revision history  (yymmdd)
c   791001  date written
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  mpunfl
      integer x(*)
c***first executable statement  mpunfl
      call mpchk (1, 4)
c the underflowing number is set to zero
c an alternative would be to call mpminr (x) and return,
c possibly updating a counter and terminating execution
c after a preset number of underflows.  action could easily
c be determined by a flag in labelled common.
      x(1) = 0
      return
      end
