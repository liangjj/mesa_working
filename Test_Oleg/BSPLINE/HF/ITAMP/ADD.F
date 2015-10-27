*
*     ----------------------------------------------------------------
*               a d d
*     ----------------------------------------------------------------
*
*     add a slater integral to the data structure associated with the
*     energy expression
*
      SUBROUTINE add(c,k,i,j,first)
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
      LOGICAL first
*
      ip = ijptr(i-nclosd,j-nclosd)
     
      if (first) then
         coef(ip+k/2+1) = c/sum(i) + coef(ip+k/2+1)
        else
           ip = ip + min(l(i),l(j)) +1 + (k-abs(l(i)-l(j)))/2 + 1
         coef(ip) = coef(ip) + c/sum(i)
        end if
      end
