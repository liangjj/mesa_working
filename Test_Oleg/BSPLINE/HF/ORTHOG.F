*
*     ------------------------------------------------------------------
*               o r t h o g
*     ------------------------------------------------------------------
*
*       this routine orthogonalizes the set of radial functions when an
*   orthogonality constraint applies.  a gram-schmidt type of  process
*   is used.
*
      SUBROUTINE orthog
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
*
      COMMON /SLATER/fc1(nt,kx),fc2(nt,kx),yk(nt),fyk(nt,kx),
     :        dyk(ktx,nt),ipvtd(nt),fc(nt,kx),ic1,ic2,iy1,iy2,ky
*
      LOGICAL fail,omit,rel,all,trace,change
      COMMON /TEST/fail,omit,rel,all,trace
*
      if (nwf .eq. 1 .or. ib .gt. nwf) return
      write (6,26)
26    format(/)
      ii = max0(2,ib)
      do 2 i = ii,nwf
         change = .false.
         do 60 j = 1,i-1
         if (e(i,j) .ne. 0.d0 ) then
*
*        orthogonality condition applies
*
            c = quadr(i,j,0)
            if (dabs(c) .gt. 1.d-10) then
                write(6,63) el(j),el(i),c
63              format(6x,'<',a3,'|',a3,'>=',1pd8.1)
                m = max0(m,max(j))
                do 64 jj = 1,m
                   p(jj,i) = p(jj,i) - c*p(jj,j)
64              continue
                change = .true.
            end if
         end if
60       continue
         if (change) then
            pnn = dsqrt(quadr(i,i,0))
         if (p(1,i) .lt. 0.d0) pnn = - pnn
         do 66 jj = 1,m
            p(jj,i) = p(jj,i)/pnn
66       continue
         m = no
67       if (dabs(p(m,i)) .lt. 1.d-15) then
            p(m,i) = 0.d0
            m = m-1
            go to 67
         end if
         max(i) = m
      end if
2     continue
      ic1 = 0
      ic2 = 0
      iy1 = 0
      iy2 = 0
      end
