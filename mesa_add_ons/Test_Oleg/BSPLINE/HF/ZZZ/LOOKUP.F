*
*     -----------------------------------------------------------------
*                l o o k - u p
*     -----------------------------------------------------------------
*
      SUBROUTINE lookup(tab,p1,p2,ind,no,key)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      INTEGER tab(*),p1,p2,ind,no,key
      do 40 i = p1,p2
         if (tab(i).eq.key) then
            no = no + 1
            ind = i
         end if
40    continue
      end
