*
*     ------------------------------------------------------------------
*               c a
*     ------------------------------------------------------------------
*
*
      double precision FUNCTION ca(l,k)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      COMMON /eav/cca(10),ccb(35)
*
      if (l .le. 4) then
         ca = cca((l*(l-1) + k)/2)
       else
         ca = rme(l,l,k)**2
      end if
      end
