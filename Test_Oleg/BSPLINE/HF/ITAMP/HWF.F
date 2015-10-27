*
*     ------------------------------------------------------------------
*               h w f
*     ------------------------------------------------------------------
*
*       returns the value of an unnormalized (nl) hydrogenic function
*   with nuclear charge zz and radius r.
*
*
      double precision FUNCTION hwf(n,l,zz,r)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
      k = n-l-1
      p = 1.d0
      a = 1.d0
      b = k
      c = n+ l
      x = -2.d0*zz*r/n
*
*  *****  TEST if underflow may occur, if so set hwf = 0
*
      if ( x .gt. -150.d0 ) then
         if (k .ge. 0) then
           do 4 i = 1,k
             p = 1.d0 + a/b*p/c*x
             a = a + 1.d0
             b = b - 1.d0
             c = c - 1.d0
4          continue
           hwf = p*dexp(x/2.d0)*(-x)**(l+1)
         else
           write(iscw,7) n,l,zz,r
7          format('forbidden combination of n and l in HWF subprogram'/
     :    'n = ',i4,'  l = ',i4,'   z = ',f6.1,'   r = ',f8.4)
           stop
	 end if
      else
        hwf = 0.d0
      end if
      end
