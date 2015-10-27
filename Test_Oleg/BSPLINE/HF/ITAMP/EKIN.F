*     ------------------------------------------------------------------
*                       e k i n
*     ------------------------------------------------------------------
*
*       returns the value of the integral of
*
*         (2/r)p (y p  + x )
*               j  i i    i
*
*   integrated with respect to r.
*
*
      double precision function ekin(i,ii,rel)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(nt=100,kx=8,nwd=10)
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      DOUBLE PRECISION yrm(nt,kx),xv(nt,kx),w(nt,kx)
      LOGICAL rel
      call xch(i,rel,xv)
      call potl(i,rel,yrm)
      call bxvpw(nt,ks,ns,2.0d0,yrm,p(1,i),xv,w)
      ekin = ddot(nt,p(1,ii),1,xv,1)
      print *, 'Ekin', i,ii,ekin
      return
      end
