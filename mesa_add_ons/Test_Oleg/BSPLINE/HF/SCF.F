*
*     ------------------------------------------------------------------
*               s c f
*     -----------------------------------------------------------------
*
*       this routine controls the scf procedure described in chapter
*   7.  if certain input parameters are zero (or blank) they will  be
*   set to their default value.
*
*          parameter       default value
*          --------        -------------
*          scftol          1.d-7
*          i*              (nwf + 1 - ib)/4 + 3
*          nscf            12
*
*   the self-consistency convergence criterion is
*
*          z2 = sqrt( scftol*(z*nwf/2) )
*
*   it is increased by a factor two at the end of each iteration.
*
*
      SUBROUTINE scf(etotal,scftol,erel)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      parameter(nt=100,kx=8,nwd=10)
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
      CHARACTER config*50,el*3,atom*6,term*6,ans*1
      COMMON /LABEL/config,el(nwd),atom,term
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      LOGICAL fail,omit,rel,all,trace,last
      COMMON /TEST/fail,omit,rel,all,trace
*
*  *****  set the scf convergence parameter to an optimistic value
*
      rel = .false.
      call energy(enonr)
      tol = dsqrt(z)*1.d-10
      z2 = scftol*dsqrt(z*nwf)
      write (6,15)
15    format(//)
      write (6,16) omit,scftol,no
   16 format(10x,44hweak orthogonalization during the scf cycle=,l4/
     :       10x,44hscf convergence tolerance (functions)      =,1pd9.2
     :      /10x,44hnumber of points in the maximum range      =,i4)
*
*  *****  set iteration parameters
*
      ipr = 0
      dp1 = 0.d0
      etotal = 0.d0
      icycle = 0
      if ( ib .gt. nwf ) go to 17
*
*  *****  perform nscf self-consistent field iterations
*
      last = .false.
9     do 100 i = 1,nscf
      icycle = icycle + 1
      write(6,7) icycle,z2
7     format(//10x,17hiteration number ,i2/10x,16h----------------//
     : 10x,'scf convergence criteria (scftol*sqrt(z*nwf)) = ',1pd9.1/)
      dp1 = 0.d0
      call grange
*
*  *****  solve each differential equation in turn
*
      write(6,14)
14    format(/20x,' el',9x,'ed',13x,'az',11x,'norm',7x,'dpm')
      do 2 jp = ib,nwf
      j = iord(jp)
      call de(j)
      if ( fail ) return
      dp = dpm(j)*dsqrt(sum(j))
      if ( dp1 .ge. dp ) go to 2
      dp1 = dp
      jj = j
2     continue
      if ( dp1 .lt. z2) go to 6
*
*  *****  solve ic differential equations each time selecting the
*  *****  one with the largest dpm
*
      do 4 ii =1,ic
      call de(jj)
      if ( fail ) return
      dp1 = 0.d0
      do 5 jp = ib,nwf
      j = iord(jp)
      dp = dsqrt(sum(j))*dpm(j)
      if ( dp1 .gt. dp ) go to 5
      jj = j
      dp1 = dp
5     continue
      if (dp1 .lt. z2) go to 6
4     continue
6     call orthog
      if (dp1 .lt. z2 .and. last) go to 17
      if ( i .eq. nscf ) go to 1
*
*  *****  if functions appear to have converged,solve each again, in
*  *****  turn, and TEST again
*
      if (dp1 .le. z2) last = .true.
*
*  *****  increase the convergence criterion for self-consistency
*
1     z2 = 2.d0*z2
      write(3,8) el(jj),dp1
      write(iscw,8) el(jj),dp1
8     format(/ 6x,34hleast self-consistent function is ,a3,
     :   27h :weighted maximum change =,1pd10.2)
100   cfgtol = 1.4d0*cfgtol
18    write(iscw,13)
13    format(10x/' scf iterations have converged to the above accuracy')
      write(3,13)
20    write(iscw,'(a)')' do you wish to continue ? (y/n/h) '
      read(5,'(a)') ans
      if (ans .eq. 'h' .or. ans .eq. 'h') then
         call help(10)
         go to 20
      end if
      if (ans .eq. 'y' .or. ans .eq. 'y') then
         write(iscw,'(a,a)')
     :     ' enter the additional iterations and new ic',
     :     ' in format(i2, 1x, i2)'
         read(5,'(i2, 1x, i2)') nscf,ic
         go to 9
      end if
      fail = .true.
*
*  *****  perform relativistic and non-relativistic calculations
*
17    continue
      call energy(enonr)
      rel = .true.
      call energy(etotal)
      rel = .false.
      erel = etotal - enonr
      nit = nwf - ib + 1
      write (3, 105) nit, dp1
105   format(//10x,'number of functions iterated          =',i6/
     :         10x,'maximum weighted change in functions  =',d10.2/)
      return
      end
