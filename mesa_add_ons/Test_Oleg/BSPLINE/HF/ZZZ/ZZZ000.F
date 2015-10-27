*
*     ------------------------------------------------------------------
*                       d e
*     ------------------------------------------------------------------
*
*       this routine controls the solution of the differenttial equation
*   for the radial function p  .  one of three methods is selected -
*                            i1
*   m1, m2, or m3 -  for solving the equations,  the  initial  choice
*   being determined by an input paramter, meth(i), except when no
*   exchange is present, in which case m2 is selected. (for further
*   information see sec. 7-4)
*
*        value of meth(i)     method
*        ---------------      ------
*        < or =1            m1 with search for an acceptable solution
*             =2            m2 with search for an acceptable solution
*             =3            m3 without any checking
*
*   if m1 fails to find an acceptable solution, the radial  functions
*   are  orthogonalized,  off-diagonal  energy parameters recomputed,
*   and the method tried again.   should it continue to fail, meth(i)
*   is set to 2.
*
*
      SUBROUTINE de(i1)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(nt=100,kx=8,ktx=3*kx-2,nwd=10)
      COMMON /coeff/coef(100),ijptr(5,5)
      CHARACTER config*50,el*3,atom*6,term*6
      COMMON /label/config,el(nwd),atom,term
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /SCF/ sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :             iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      COMMON /SLATER/fc1(nt,kx),fc2(nt,kx),yk(nt),fyk(nt,kx),
     :        dyk(ktx,nt),ipvtd(nt),fc(nt,kx),ic1,ic2,iy1,iy2,ky
*
      LOGICAL fail,omit,rel,all,trace,change
      COMMON /test/fail,omit,rel,all,trace
      DOUBLE PRECISION pde(nt)
      CHARACTER*3 aster(3)
      LOGICAL first
      data aster(1),aster(2),aster(3)/'  ','* ','**'/
      DATA ipr,ed1,ed2/0,0.d0,0.d0/
*
      i = i1
      ed2 = e(i,i)
      mk= max0(1,meth(i))
      if (nwf .eq. 1) mk = 2
      node = n(i) - l(i) - 1
      first = .true.
*
*  *****  call method to solve the differential equation
*
12      call method(i,first,pde,m)
        if ( fail ) then
*
*         .. mk failed to find an acceptable solution, orthogonalize
*            the estimates and try again if i ge. ib
          if (i .ge. ib .and. first) then
            call orthog
            call grange
	    first = .false.
            call method(i,first,pde,m)
	  end if
          if ( fail ) then
*
*           .. error return from second try.  if m1 was used,switch to
*              m2 and try again
            if ( mk .eq. 2) return
            mk = 2
            go to 12
	  end if
        end if
        pn = dsqrt(quad(ns,pde,pde))
        do 9 j = 1,ns
          pde(j) = pde(j)/pn
9       continue
*
*  *****  check if different method should be used
*
        if ( mk .eq. 1 ) then
          if (dabs(1.d0 -ed2/e(i,i)) .lt. 0.005d0  .and.
     :       dmax1(dabs(1.d0-pn), dabs(1.d0/pn-1.d0)).gt. .20d0 ) then
            meth(i) = 2
            mk = 2
            go to 12
          end if
        else
          if (dabs(1.d0 - ed2/e(i,i)) .lt. 0.0001d0 .and.
     :       ic .gt. 1) ic = ic -1
        end if
*
*  *****  set the accelerating parameter
*
13      if (ipr .ne. i) then
          acc(i) = .75*acc(i)
        else
          ed2 = ed2 - e(i,i)
          if (ed1*ed2 .gt. 0.d0 ) then
            acc(i) = .75*acc(i)
          else
            acc(i) = (1.d0 + 3.d0*acc(i))/4.d0
          end if
        end if
        c = acc(i)
        cd = 1.d0 - c
*
*   *****  improve the estimates
*
        max(i) = m
        dp     = 0.d0
        do 21 j = 1,m
          diff = p(j,i)-pde(j)
          dp     = dmax1(dp    ,dabs(diff))
          p(j,i) = pde(j) + c*diff
21     continue
       if (m .lt. ns) then
         m = m + 1
         do 24 j = m,ns
           p(j,i) = 0.d0
24       continue
       end if
*
*  *****  check the orthogonalization
*
26    nn = nwf
      if (omit) nn = ib - 1
      ij = 0
      dpw = dp/dsqrt(sum(i))
      m = max(i)
      change = .false.
      do 60 j = 1,nn
         if (e(i,j) .ne. 0.d0 .and. i .ne. j ) then
          if (dpm(j) .lt. dsqrt(sum(j))*dpw .or. j .lt. ib ) then
*
*        orthogonality condition applies
*
            c = quadr(i,j,0)
            write(6,63) el(j),el(i),c
63          format(6x,'<',a3,'|',a3,'>=',1pd8.1)
            m = max0(m,max(j))
            do 64 jj = 1,m
               p(jj,i) = p(jj,i) - c*p(jj,j)
64          continue
            change = .true.
          end if
       end if
60    continue
      if (change .or. c.ne.0.d0) then
         pnn = dsqrt(quadr(i,i,0))
         do 66 j = 1,m
            p(j,i) = p(j,i)/pnn
66       continue
      end if
67    if (dabs(p(m,i)) .lt. 1.d-15) then
         p(m,i) = 0.d0
         m = m-1
         go to 67
      end if
      max(i) = m
      write (6,17) el(i),e(i,i),az(i),pn,aster(mk),dp
17    format(20x,a3,2f15.7,f12.7, a2,1pd10.2)
      dpm(i) = dp
      if (ipr .eq. i1) then
         ed1 = ed2
      else
         ed1 = ed2 - e(i1,i1)
      end if
*     .. reset the Slater integral parameters
      ipr = i1
      ic1 = 0
      ic2 = 0
      iy1 = 0
      iy2 = 0
      end
