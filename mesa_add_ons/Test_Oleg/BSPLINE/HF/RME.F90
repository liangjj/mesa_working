!=======================================================================
  FUNCTION rme(l,lp,k) 
!=======================================================================
!  Evaluates the reduced matrix element (l//c(k)//lp)  -  see fano
!    and racah, irreducible tensorial sets, chap. 14, p. 81
!
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l,lp,k
    REAL(KIND=8) :: rme, qusqrt
    INTEGER :: i2g, ig, i1,i2,i3

      if (min(l,lp) == 0) then
         rme = 1
       else if ( k == 0) then
         rme = 2*l+1
         rme = dsqrt(rme)
       else
         i2g=l+lp+k
         ig=i2g/2
         if (i2g - 2*ig <> 0) then
             rme = 0
         else
            i1=ig-l
            i2=ig-lp
            i3=ig-k
            qusqrt=(2*l+1)*(2*lp+1)
            rme=dsqrt(qusqrt)*dexp((gam(2*i1+1)+gam(2*i2+1)+gam(2*i3+1)-
     :        gam(i2g+2))/2.d0 +gam(ig+1)-gam(i1+1)-gam(i2+1)-gam(i3+1))
         end if
      end if
   END FUNCTION rme
