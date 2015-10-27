*deck mpmlp
      subroutine mpmlp (u, v, w, j)
c***begin prologue  mpmlp
c***subsidiary
c***purpose  subsidiary to dqdota and dqdoti
c***library   slatec
c***type      all (mpmlp-a)
c***author  (unknown)
c***description
c
c performs inner multiplication loop for mpmul. carries are not pro-
c pagated in inner loop, which saves time at the expense of space.
c
c***see also  dqdota, dqdoti
c***routines called  (none)
c***revision history  (yymmdd)
c   791001  date written
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  mpmlp
      integer u(*), v(*), w
c***first executable statement  mpmlp
      do 10 i = 1, j
   10 u(i) = u(i) + w*v(i)
      return
      end
