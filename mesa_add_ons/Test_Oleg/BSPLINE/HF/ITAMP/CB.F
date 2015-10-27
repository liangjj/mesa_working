*
*     -----------------------------------------------------------------
*                 c b
*     -----------------------------------------------------------------
*
*
      double precision FUNCTION cb(l,lp,k)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      COMMON /eav/cca(10),ccb(35)
      INTEGER icbptr(0:4)
      data    icbptr/1,6,14,23,31/
*
             if (l .le. lp) then
                 l1 = l
                 l2 = lp
              else
                 l1 = lp
                 l2 = l
             end if
             if ( l2 .le. 4) then
                cb = ccb(icbptr(l1)+(k+l1-l2)/2+(l1+1)*(l2-l1))
               else
                cb = rme(l,lp,k)**2/(2*(2*l+1)*(2*lp+1))
             end if
      end
