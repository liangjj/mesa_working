*
*
*     ------------------------------------------------------------------
*               h n o r m
*     ------------------------------------------------------------------
*
*       returns the value of the normalization constant for an (nl)
*   hydrogenic function with nuclear charge zz.
*
*
      double precision FUNCTION hnorm(n,l,zz)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      m = l + l + 1
      a = n + l
      b = m
      t = a
      d = b
      m = m - 1
      do 1 i = 1,m
        a = a - 1.d0
        b = b - 1.d0
        t = t*a
        d = d*b
   1  continue
      hnorm = dsqrt(zz*t)/( n*d)
      end
