*deck crotg
      subroutine crotg (ca, cb, c, s)
c***begin prologue  crotg
c***purpose  construct a givens transformation.
c***library   slatec (blas)
c***category  d1b10
c***type      complex (srotg-s, drotg-d, crotg-c)
c***keywords  blas, givens rotation, givens transformation,
c             linear algebra, vector
c***author  (unknown)
c***description
c
c    complex givens transformation
c
c    construct the givens transformation
c
c             (c    s)
c       g  =  (      ),  c**2 + abs(s)**2 =1,
c             (-s   c)
c
c    which zeros the second entry of the complex 2-vector (ca,cb)**t
c
c    the quantity ca/abs(ca)*norm(ca,cb) overwrites ca in storage.
c
c    input:
c        ca (complex)
c        cb (complex)
c
c    output:
c        ca (complex)      ca/abs(ca)*norm(ca,cb)
c        c  (real)
c        s  (complex)
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   790101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  crotg
      complex ca, cb, s
      real c
      real norm, scale
      complex alpha
c***first executable statement  crotg
      if (abs(ca) .eq. 0.0) then
        c = 0.0
        s = (1.0,0.0)
        ca = cb
      else
        scale = abs(ca) + abs(cb)
        norm = scale * sqrt((abs(ca/scale))**2 + (abs(cb/scale))**2)
        alpha = ca /abs(ca)
        c = abs(ca) / norm
        s = alpha * conjg(cb) / norm
        ca = alpha * norm
      endif
      return
      end
