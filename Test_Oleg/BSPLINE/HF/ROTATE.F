*
*     ------------------------------------------------------------------
*               r o t a t e
*     ------------------------------------------------------------------
*
*        this routine analyses the energy expression to determine the
*   stationary condition with respect to rotation of orbials i and j.
*   if the condition is zero, the off-diagonal energy parameters may
*   be set to zero;  otherwise the orbials are rotated so as to satisfy
*   the stationay condition to first order in the perturbation.
*
*
      SUBROUTINE rotate(i,j)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(nt=100,kx=8,ktx=3*kx-2,nwd=10)
      COMMON /coeff/coef(100),ijptr(5,5)
      CHARACTER config*50,el*3,atom*6,term*6
      COMMON /LABEL/config,el(nwd),atom,term
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff

      COMMON /SLATER/fc1(nt,kx),fc2(nt,kx),yk(nt),fyk(nt,kx),
     :        dyk(ktx,nt),ipvtd(nt),fc(nt,kx),ic1,ic2,iy1,iy2,ky
*
      LOGICAL fail,omit,rel,all,trace
      COMMON /TEST/fail,omit,rel,all,trace
      all = .true.
      g = 0.d0
      dg = 0.d0
      qi = sum(i)
      qj = sum(j)
      if (qi .eq. 2.d0*(2*l(i)+1) .and. qj .eq. 2.d0*(2*l(j)+1))
     :   go to 44
      if (dabs(qi - qj) .gt. 1.d-14) then
        c = 0.5d0*(qi - qj)
        g = g -c*hl(el,i,j,rel)
        dg = dg -c*(hl(el,i,i,rel) - hl(el,j,j,rel))
      end if
*
16    do 13 k = 0,2*l(i),2
        c = qi*(a(i,i,k) - a(i,j,k) - b(i,j,k))
        if (dabs(c) .gt. 1.d-8) then
          g = g  + c*rky(i,i,i,j,k,rel)
          fkii = fky(i,i,k,rel)
          fkij = fky(i,j,k,rel)
          gkij = gky(i,j,k,rel)
          dg = dg +c*(fkii - fkij - 2.d0*gkij)
        end if
        cj = qj*(a(j,j,k) - a(j,i,k) - b(j,i,k))
        if (dabs(cj) .gt. 1.d-8) then
          fkjj = fky(j,j,k,rel)
          if ( .not.(dabs(c) .gt. 1.d-8)) then
            fkij = fky(i,j,k,rel)
            gkij = gky(i,j,k,rel)
	end   if
          g = g - cj*rky(j,j,j,i,k,rel)
          dg = dg + cj*(fkjj -fkij - 2.d0*gkij)
        end if
13    continue
      do 14 m = 1,nwf
        if (m .eq. i .or. m.eq. j) go to 14
        do 15 k = 0,2*min0(l(i),l(m)),2
          c = a(i,m,k)*qi - a(j,m,k)*qj
          if (dabs(c) .gt. 1.d-8) then
            g = g + c*rky(i,m,j,m,k,rel)
            dg = dg + c*(fky(i,m,k,rel) - fky(j,m,k,rel))
	  end if
23        kk = iabs(l(i)-l(m)) + k
          c = b(i,m,kk)*qi - b(j,m,kk)*qj
          if  (dabs(c) .gt. 1.d-8) then
            g = g + c*rky(i,j,m,m,kk,rel)
            dg = dg + c*(gky(i,m,kk,rel) - gky(j,m,kk,rel))
	  end if
15      continue
14    continue
      if (dabs(qi-qj) + dabs(g) + dabs(dg).le. 2.d-8 ) go to 44
      print *, 'qi-qj',qi-qj
      print *, 'g,dg', g,dg
      print *, 'e(i,j)',i,j,e(i,j)
101   if (dabs(g)+dabs(dg) .gt. 1.d-8 .or. dabs(e(i,j)) .gt. 2.d-5) then
         eps = g/dg
         eps = dsign(dmin1(dabs(eps),0.2d0),eps)
         dd = dsqrt(1.d0 + eps*eps)
	 print *, 'eps,dd,ns',eps,dd,ns
         do 41 jj = 1,ns
            pi = (p(jj,i) + eps*p(jj,j))/dd
            p(jj,j) = (p(jj,j) - eps*p(jj,i))/dd
41          p(jj,i) = pi
	 ic1 = 0 
	 ic2 = 0
	 iy1 = 0
	 iy2 = 0
      else
         eps = 0.d0
      end if
      write (6,100) el(i),el(j),g,el(i),el(j),dg,eps
100   format(10x,'c(',2a3,') =',f12.5,3x,'v(',2a3,') =',f12.5,
     :       3x,'eps =',f9.6)
      return
*
*  *****  the energy is stationary with respect to rotations
*
44    e(i,j) = 1.d-10
      e(j,i) = 1.d-10
      return
      end
