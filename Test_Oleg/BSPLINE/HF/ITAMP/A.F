*
*     ----------------------------------------------------------------
*               a
*     ----------------------------------------------------------------
*
*       determine the coefficient in the potential for electron i of
*       y^k(j,j)
     
      double precision FUNCTION a(i,j,k)
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
      if (i.gt.nclosd .and. j.gt.nclosd) then
         istart = ijptr(i-nclosd,j-nclosd) + 1
         a = coef(istart + k/2)
      else if (i.eq.j) then
         c = sum(i) - 1.d0
         if (k.eq.0) then
            a = c
         else
            a = -c*ca(l(i),k)
         end if
      else if (k.eq.0) then
         a = sum(j)
      else
         a = 0.d0
      end if
      end
