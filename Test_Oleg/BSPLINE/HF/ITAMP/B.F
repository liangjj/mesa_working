*
*     ----------------------------------------------------------------
*               b
*     ----------------------------------------------------------------
*
*     determine the coefficient of the y^k(i,j)p(j) term in the exchange
*     expression of electron i
*
      double precision FUNCTION b(i,j,k)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(nt=100,kx=8,nwd=10)
      COMMON /coeff/coef(100),ijptr(5,5)
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      if (i.eq.j) then
         b = 0.d0
      else if (i.gt.nclosd .and. j.gt.nclosd) then
*
*   ..... ll is the number of direct terms
*         istart the beginning of the exchange terms
*
         ll = min(l(i),l(j)) + 1
         istart = ijptr(i-nclosd,j-nclosd) + 1 + ll
         kk = (k - abs(l(i)-l(j)))/2
         b = coef(istart + kk)
      else
         b = -sum(j)*cb(l(i),l(j),k)
      end if
      end
