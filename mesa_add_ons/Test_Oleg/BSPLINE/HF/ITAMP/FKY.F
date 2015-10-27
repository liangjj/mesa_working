*
*     ------------------------------------------------------------------
*              F K Y
*     ------------------------------------------------------------------
*                             k
*       Returns the value of F (i,j)
*
      DOUBLE PRECISION FUNCTION FKy(I,J,K,rel)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      fky = rky(i,j,i,j,k,rel)
      end
