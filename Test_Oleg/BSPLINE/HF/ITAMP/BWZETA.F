*
*     ----------------------------------------------------------------
*               b w z e t a
*     ----------------------------------------------------------------
*
*
*  ***** computes the nuclear spin-orbit parameter and the
*        corrections for the other electrons
*      using the formula derived by blume and watson.
*
      double precision FUNCTION bwzeta(i1)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(nt=100,kx=8,nwd=10)
      COMMON/blume/coefn2(4),coefnk(4),coefvk(4)
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      dimension ss(3)
*
*
      zeta = fine*z*quadr(i1,i1,-3)
      lb = l(i1)
      do 10 i = 1,nwf
       if (i .eq. i1) go to 10
         la = l(i)
         zeta = zeta -sum(i)*sn(i1, i, i1, i, 0)
       if (sum(i) .ne. 4*l(i)+2) go to 10
         call bwint(la,lb)
         ke1 = 2
         if (la .ne. lb) ke1 = iabs(la-lb)
         ip = 0
         do 20 k = ke1,la+lb,2
            ip = ip+1
            zeta = zeta+coefn2(ip)*sn(i1, i, i, i1, k-2)
     :                 +coefnk(ip)*sn(i, i1, i1, i, k)
     :                 +coefvk(ip)*(vk(i1,i,i,i1,k-1)-vk(i,i1,i1,i,k-1))
   20    continue
   10 continue
      zeta = 2.d0*zeta
      c= sum(i1)
      if (c .ne. 1.d0) then
         ss(1) = sn(i1,i1,i1,i1,0)
         c = c + c - 3.d0
         zeta = zeta - c*ss(1)
         if (lb .eq. 2) then
            ss(2) = sn(i1,i1,i1,i1,2)
            zeta = zeta + ss(2)*6.d0/7.d0
         else if (lb .eq. 3) then
            ss(2) = sn(i1,i1,i1,i1,2)
            ss(3) = sn(i1,i1,i1,i1,4)
            zeta = zeta + ss(2) + ss(3)/2.2d0
         end if
      end if
      bwzeta = zeta
      end
