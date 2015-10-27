*     ==================================================================
*              G K
*     ==================================================================
*                             k
*       Returns the value of G (i,j).
*
*
      DOUBLE PRECISION FUNCTION GKy(I,J,K,REL)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      gky = rky(i,j,j,i,k,rel)
      end
