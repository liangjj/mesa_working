*
*     -----------------------------------------------------------------
*           f a c t r l
*     -----------------------------------------------------------------
*
*
      SUBROUTINE factrl(nfact)
*
*      gam(i) = log( gamma(i-1) ), where gamma(i) = factorial i-1
*
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
*
      COMMON/fact/gam(100)
      data zero,one,two/0.d0,1.d0,2.d0/
*
      gamma=one
      gam(1) = zero
      do 1 i=1,nfact-1
         gamma=i*gamma
         gam(i+1) = dlog(gamma)
    1 continue
      do 20 i = nfact+1,(100)
         x = i-1
         gam(i) = gam(i-1) + dlog(x)
   20 continue
      return
      end
