!=======================================================================
   SUBROUTINE de(i1) 
!=======================================================================
!
!
!       This routine controls the solution of the differenttial equation
!   for the radial function p  .  one of three methods is selected -
!                            i1
!   m1, m2, or m3 -  for solving the equations,  the  initial  choice
!   being determined by an input paramter, meth(i), except when no
!   exchange is present, in which case m2 is selected. (for further
!   information see sec. 7-4)
!
!        value of meth(i)     method
!        ---------------      ------
!        < or =1            m1 with search for an acceptable solution
!             =2            m2 with search for an acceptable solution
!             =3            m3 without any checking
!
!   if m1 fails to find an acceptable solution, the radial  functions
!   are  orthogonalized,  off-diagonal  energy parameters recomputed,
!   and the method tried again.   should it continue to fail, meth(i)
!   is set to 2.
!
!----------------------------------------------------------------------
!
    IMPLICIT NONE
      PARAMETER(nt=100,kx=8,ktx=3*kx-2,nwd=10)
      COMMON /coeff/coef(100),ijptr(5,5)
      CHARACTER config*50,el*3,atom*6,term*6
      COMMON /label/config,el(nwd),atom,term
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/ sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :             iord(nwd), e(nwd,nwd)
!
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
!
      COMMON /SLATER/fc1(nt,kx),fc2(nt,kx),yk(nt),fyk(nt,kx),
     :        dyk(ktx,nt),ipvtd(nt),fc(nt,kx),ic1,ic2,iy1,iy2,ky
!
      LOGICAL fail,omit,rel,all,trace,change
      COMMON /test/fail,omit,rel,all,trace
      DOUBLE PRECISION pde(nt)
      CHARACTER*3 aster(3)
      LOGICAL first
      data aster(1),aster(2),aster(3)/'  ','* ','**'/
      DATA ipr,ed1,ed2/0,0.d0,0.d0/
!
      i = i1
      ed2 = e(i,i)
      mk= max0(1,meth(i))
      if (nwf == 1) mk = 2
      first = .true.
!
!  *****  call method to solve the differential equation
!
12      call method(i,first,ed,pde,ms)
        if ( fail ) then
!
!         .. mk failed to find an acceptable solution, orthogonalize
!            the estimates and try again if i ge. ib
          if (i .ge. ib .and. first) then
            call orthog
            call grange
	    first = .false.
            call method(i,first,ed,pde,ms)
	  end if
          if ( fail ) then
!
!           .. error return from second try.  if m1 was used,switch to
!              m2 and try again
            if ( mk .eq. 2) return
            mk = 2
            go to 12
	  end if
        end if
        pn = dsqrt(quad(ns,ks,pde,pde))
        pde = pde/pn
	e(i,i) = ed
!
!  *****  check if different method should be used
!
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
!
!  *****  set the accelerating parameter
!
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
!
!   *****  improve the estimates
!
        max(i) = ms
        dp     = 0.d0
	print *, 'P(.i)'
	print '(6f12.8)', (p(ii,i),ii=1,ns)
	print *, 'pde'
	print '(6f12.8)', (pde(ii),ii=1,ns)
        do 21 j = 1,ms
          diff = p(j,i)-pde(j)
          dp     = dmax1(dp    ,dabs(diff))
          p(j,i) = pde(j) + c*diff
21     continue
       if (ms .lt. ns) then
         ms = ms + 1
         p(ms:ns,i) = 0.d0
       end if
!
!  *****  check the orthogonalization
!
26    nn = nwf
!     if (omit) nn = ib - 1
      ij = 0
      dpw = dp/dsqrt(sum(i))
      ms = max(i)
      change = .false.
      do 60 j = 1,nn
!        if (e(i,j).ne.0.d0 .and. e(i,j).ne. 1.d-10 .and. i.ne.j ) then
         if (e(i,j).ne.0.d0 .and.  i.ne.j ) then
          if (dpm(j) .lt. dsqrt(sum(j))*dpw .or. j .lt. ib ) then
!
!        orthogonality condition applies
!
            c = quadr(i,j,0)
            write(6,63) el(j),el(i),c
63          format(6x,'<',a3,'|',a3,'>=',1pd8.3)
            ms = max0(ms,max(j))
            p(1:ms,i) = p(1:ms,i) - c*p(1:ms,j)
	    print *, 'Check orthogonality',quadr(i,j,0)
            change = .true.
	  else
!           .. the other function should be orthogonalized
            c = quadr(i,j,0)
	    pn = 1.d0/sqrt(1.d0 - c*c)
            write(6,63) el(i),el(j),c
            ms = max0(ms,max(j))
            p(1:ms,j) = (p(1:ms,j) - c*p(1:ms,i))*pn
	    print *, 'Check orthogonality',quadr(i,j,0)
          end if
       end if
60    continue
      if (change ) then
         pnn = 1.d0/sqrt(quadr(i,i,0))
         p(1:ms,i) = p(1:ms,i)*pnn
      end if
67    if (dabs(p(ms,i)) .lt. 1.d-15) then
         p(ms,i) = 0.d0
         ms = ms-1
         go to 67
      end if
      max(i) = ms+1
      azz = az(z,h,ks,l(i)+1,p(l(i)+2,i))
      write (6,17) el(i),e(i,i),azz,pn,aster(mk),dp
17    format(20x,a3,2f15.7,f12.7, a2,1pd10.2)
      dpm(i) = dp
      if (ipr .eq. i1) then
         ed1 = ed2
      else
         ed1 = ed2 - e(i1,i1)
      end if
      print *, 'ed1,ed2',ed1,ed2
!     .. reset the Slater integral parameters
      ipr = i1
      ic1 = 0
      ic2 = 0
      iy1 = 0
      iy2 = 0
      print *, ' New orbital',i,quadr(1,2,0)
      print '(6f12.8)', (p(jj,i),jj=1,ns)
      END SUBROUTINE de
