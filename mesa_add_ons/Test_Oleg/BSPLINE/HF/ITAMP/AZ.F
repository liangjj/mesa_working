*     ==================================================================
*          A Z
*     ==================================================================
*
*      Value of B_lp/r^lp at r = 0 where lp = l+1, and p = coefficient
*      of B_lp(r) in the expansion of the orbital.
*
       DOUBLE PRECISION FUNCTION az(z,h,ks,lp,p)
       IMPLICIT DOUBLE PRECISION (a-h,o-z)
*
	if (lp .lt. ks ) then
	   coeff = 1.d0
	   c = z/h
	   do 20 j = 1,lp
	     coeff = (coeff*c*(ks-j))/(j*j)
   20     continue
	else
	  coeff = 0.d0
	end if
   10   continue
	az = coeff*p
	end
