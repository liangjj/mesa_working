*     ------------------------------------------------------------------
*                r m e
*     ------------------------------------------------------------------
*
*
      double precision FUNCTION rme(l,lp,k)
*
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
*
      COMMON/fact/gam(100)
*
*--- evaluates the reduced matrix element (l//c(k)//lp)  -  see fano
*    and racah, irreducible tensorial sets, chap. 14, p. 81
*
      if (min0(l,lp) .eq. 0) then
         rme = 1.d0
       else if ( k .eq. 0) then
         rme = 2*l+1
         rme = dsqrt(rme)
       else
         i2g=l+lp+k
         ig=i2g/2
         if (i2g - 2*ig .ne. 0) then
             rme = 0.d0
         else
            i1=ig-l
            i2=ig-lp
            i3=ig-k
            qusqrt=(2*l+1)*(2*lp+1)
            rme=dsqrt(qusqrt)*dexp((gam(2*i1+1)+gam(2*i2+1)+gam(2*i3+1)-
     :        gam(i2g+2))/2.d0 +gam(ig+1)-gam(i1+1)-gam(i2+1)-gam(i3+1))
         end if
      end if
      return
      end
