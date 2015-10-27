*deck mpadd
      subroutine mpadd (x, y, z)
c***begin prologue  mpadd
c***subsidiary
c***purpose  subsidiary to dqdota and dqdoti
c***library   slatec
c***type      all (mpadd-a)
c***author  (unknown)
c***description
c
c adds x and y, forming result in z, where x, y and z are 'mp'
c  (multiple precision) numbers.  four guard digits are used,
c  and then r*-rounding.
c
c***see also  dqdota, dqdoti
c***routines called  mpadd2
c***revision history  (yymmdd)
c   791001  date written
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  mpadd
      integer x(*), y(*), z(*)
c***first executable statement  mpadd
      call mpadd2 (x, y, z, y, 0)
      return
      end
