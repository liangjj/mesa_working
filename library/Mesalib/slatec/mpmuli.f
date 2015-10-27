*deck mpmuli
      subroutine mpmuli (x, iy, z)
c***begin prologue  mpmuli
c***subsidiary
c***purpose  subsidiary to dqdota and dqdoti
c***library   slatec
c***type      all (mpmuli-a)
c***author  (unknown)
c***description
c
c multiplies 'mp' x by single-precision integer iy giving 'mp' z.
c this is faster than using mpmul.  result is rounded.
c multiplication by 1 may be used to normalize a number
c even if the last digit is b.
c
c***see also  dqdota, dqdoti
c***routines called  mpmul2
c***revision history  (yymmdd)
c   791001  date written
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  mpmuli
      integer x(*), z(*)
c***first executable statement  mpmuli
      call mpmul2 (x, iy, z, 0)
      return
      end
