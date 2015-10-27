!=======================================================================
  FUNCTION a(i,j,k)
!=======================================================================
!
!  determine the coefficient in the potential for electron i of
!  y^k(j,j)
!
!----------------------------------------------------------------------
!
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: i,j,k
     REAL(KIND=8) :: a

      PARAMETER(nt=100,kx=8,nwd=10)
      COMMON /coeff/coef(100),ijptr(5,5)
 
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
 
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
 
      if (i > nclosd .and. j > nclosd) then
         istart = ijptr(i-nclosd,j-nclosd) + 1
         a = coef(istart + k/2)
      else if (i == j) then
         c = sum(i) - 1.d0
         if (k == 0) then
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
!=======================================================================
  SUBROUTINE add(c,k,i,j,first) 
!=======================================================================
!   add a slater integral to the data structure associated with the
!   energy expression
!
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k,i,j
    LOGICAL, INTENT(IN) :: first
    REAL(KIND=8), INTENT(IN) :: c

    ! .. local variables
      PARAMETER(nt=100,kx=8,nwd=10)
      COMMON /coeff/coef(100),ijptr(5,5)
 
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
 
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
 
      ip = ijptr(i-nclosd,j-nclosd)
     
      if (first) then
         coef(ip+k/2+1) = c/sum(i) + coef(ip+k/2+1)
        else
           ip = ip + min(l(i),l(j)) +1 + (k-abs(l(i)-l(j)))/2 + 1
         coef(ip) = coef(ip) + c/sum(i)
        end if

  END FUNCTION add
!=======================================================================
  SUBROUTINE array 
!=======================================================================
!
!----------------------------------------------------------------------
!
    IMPLICIT NONE

      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(nt=100,kx=8,nwd=10)
      COMMON /coeff/coef(100),ijptr(5,5)
 
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
 
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
 
      ip = 0
      do i = nclosd+1,nwf
        isumi = sum(i)
        dsumi = sum(i) - isumi
        do j = nclosd+1,nwf
          isumj = sum(j)
          dsumj = sum(j) - isumj
          if ( i .ne. j) then
            c = sum(j)
            if (dsumi .ne. 0.d0 .and. dsumj .ne. 0.d0)
     :        c = (dsumi*(isumi+1)*isumj +
     :             dsumj*(isumj+1)*isumi)/sum(i)
            else
              c = sum(i) - 1.d0
              if (dsumi .ne. 0.d0)
     :          c = (isumi*(sum(i)+dsumi-1))/sum(i)
            end if
 
            ijptr(i-nclosd,j-nclosd) = ip
 
!           ...direct contribution
 
            do k = 0,2*min0(l(i),l(j)),2
              ip = ip + 1
              if (ip .gt. (100))
     :          stop ' coef array too small: max = (100)'
              coef(ip) = 0.d0
              if (k .eq. 0) then
                coef(ip) = c
              else if (i .eq. j) then
                coef(ip) = -c*ca(l(i),k)
              end if
            continue
 
!           ... exchange contribution
 
            if (i .ne. j) then
              do k = abs(l(i)-l(j)),l(i)+l(j),2
                ip = ip + 1
                if (ip .gt. (100))
     :            stop ' coef array too small: max = (100)'
                coef(ip) = -c*cb(l(i),l(j),k)
              continue
            end if
          continue
        continue

      END SUBROUTINE array
!=======================================================================
  FUNCTION b(i,j,k) 
!=======================================================================
!
!   determine the coefficient of the y^k(i,j)p(j) term in the exchange
!   expression of electron i
!
!----------------------------------------------------------------------
!
    IMPLICIT NONE

 
      PARAMETER(nt=100,kx=8,nwd=10)
      COMMON /coeff/coef(100),ijptr(5,5)
 
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
 
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
 
      if (i == j) then
         b = 0.d0
      else if (i > nclosd .and. j > nclosd) then
 
    ..... ll is the number of direct terms
          istart the beginning of the exchange terms
 
         ll = min(l(i),l(j)) + 1
         istart = ijptr(i-nclosd,j-nclosd) + 1 + ll
         kk = (k - abs(l(i)-l(j)))/2
         b = coef(istart + kk)
      else
         b = -sum(j)*cb(l(i),l(j),k)
      end if
 
  END FUNCTION b
!=======================================================================
  SUBROUTINE bwint(lc,lo)
  
!=======================================================================
!   Compute the spin-orbit parameter using the Blume and Watson method
!
!----------------------------------------------------------------------
!
    IMPLICIT NONE
      COMMON/blume/coefn2(4),coefnk(4),coefvk(4)
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
!
! ... lc is the l-value of the filled subshell, lo is the l-value
!     of the partially-filled subshell.
!
      if(lc > 3 .or. lo >  4) then
      write(iscw,'(A,I3,A,I3)') 'Incorrect calling of bwint with lc =', &`j
         lc, 'and lo =',lo
         return
      end if
      lc1 = lc + 1

      go to (10,20,30,40), lc1
   10 go to (11,12,13,14), lo
!
! ... s-p
!
   11 coefnk(1) = 1.d0
      coefn2(1) = -2.d0
      coefvk(1) = 1.d0
      return
!
! ... s-d
!
   12 coefnk(1) = 6.d0/5.d0
      coefn2(1) = -9.d0/5.d0
      coefvk(1) = 3.d0/5.d0
      return
!
! ... s-f
!
   13 coefnk(1) = 9.d0/7.d0
      coefn2(1) = -12.d0/7.d0
      coefvk(1) = 3.d0/7.d0
      return
!
! ... s-g
!
   14 coefnk(1) = 4.d0/3.d0
      coefn2(1) = -5.d0/3.d0
      coefvk(1) = 1.d0/3.d0
      return
   20 go to (21,22,23,24), lo
!
! ... p-p
!
   21 coefnk(1) = 0.d0
      coefn2(1) = 3.d0
      coefvk(1) = 9.d0/5.d0
      return
!
! ... p-d
!
   22 coefnk(1) = 3.d0/7.d0
      coefnk(2) = 36.d0/35.d0
      coefn2(1) = -12.d0/5.d0
      coefn2(2) = 0.d0
      coefvk(1) = 3.d0/5.d0
      coefvk(2) = 36.d0/35.d0
      return
!
! ... p-f
!
   23 coefnk(1) = 1.d0/7.d0
      coefnk(2) = 10.d0/7.d0
      coefn2(1) = -18.d0/7.d0
      coefn2(2) = 0.d0
      coefvk(1) = 18.d0/35.d0
      coefvk(2) = 5.d0/7.d0
      return
!
! ... p-g
!
   24 coefnk(1) = 5.d0/77.d0
      coefnk(2) = 18.d0/11.d0
      coefn2(1) = -18.d0/7.d0
      coefn2(2) = 0.d0
      coefvk(1) = 3.d0/7.d0
      coefvk(2) = 6.d0/11.d0
      return
   30 go to (31,32,33,34), lo
!
! ... d-p
!
   31 coefnk(1) = 59.d0/7.d0
      coefnk(2) = -18.d0/7.d0
      coefn2(1) = -4.d0
      coefn2(2) = 0.d0
      coefvk(1) = -1.d0
      coefvk(2) = 18.d0/7.d0
      return
!
! ... d-d
!
   32 coefnk(1) = 6.d0/7.d0
      coefnk(2) = 0.d0
      coefn2(1) = 3.d0
      coefn2(2) = 0.d0
      coefvk(1) = 3.d0/7.d0
      coefvk(2) = 10.d0/7.d0
      return
!
! ... d-f
!
   33 coefnk(1) = 9.d0/7.d0
      coefnk(2) = -13.d0/77.d0
      coefnk(3) = 75.d0/77.d0
      coefn2(1) = -18.d0/7.d0
      coefn2(2) = 0.d0
      coefn2(3) = 0.d0
      coefvk(1) = 3.d0/7.d0
      coefvk(2) = 3.d0/7.d0
      coefvk(3) = 75.d0/77.d0
      return
!
! ... d-g
!
   34 coefnk(1) = 741.d0/693.d0
      coefnk(2) = -215.d0/429.d0
      coefnk(3) = 210.d0/143.d0
      coefn2(1) = -3.d0
      coefn2(2) = 0.d0
      coefn2(3) = 0.d0
      coefvk(1) = 3.d0/7.d0
      coefvk(2) = 255.d0/693.d0
      coefvk(3) = 105.d0/143.d0
      return
   40 go to (41,42,43,44), lo
!
! ... f-p
!
   41 coefnk(1) = 52.d0/3.d0
      coefnk(2) = -20.d0/3.d0
      coefn2(1) = -9.d0
      coefn2(2) = 0.d0
      coefvk(1) = -9.d0/5.d0
      coefvk(2) = 10.d0/3.d0
      return
!
! ... f-d
!
   42 coefnk(1) = 5.d0
      coefnk(2) = 142.d0/55.d0
      coefnk(3) = -20.d0/11.d0
      coefn2(1) = -18.d0/5.d0
      coefn2(2) = 0.d0
      coefn2(3) = 0.d0
      coefvk(1) = -3.d0/5.d0
      coefvk(2) = 2.d0/5.d0
      coefvk(3) = 20.d0/11.d0
      return
!
! ... f-f
!
   43 coefnk(1) = 1.d0
      coefnk(2) = 5.d0/11.d0
      coefnk(3) = 0.d0
      coefn2(1) = 3.d0
      coefn2(2) = 0.d0
      coefn2(3) = 0.d0
      coefvk(1) = 1.d0/5.d0
      coefvk(2) = 5.d0/11.d0
      coefvk(3) = 175.d0/143.d0
      return
!
! ... f-g
!
   44 coefnk(1) = 53.d0/33.d0
      coefnk(2) = 57.d0/143.d0
      coefnk(3) = -115.d0/429.d0
      coefnk(4) = 392.d0/429.d0
      coefn2(1) = -8.d0/3.d0
      coefn2(2) = 0.d0
      coefn2(3) = 0.d0
      coefn2(4) = 0.d0
      coefvk(1) = 1.d0/3.d0
      coefvk(2) = 3.d0/11.d0
      coefvk(3) = 57.d0/143.d0
      coefvk(4) = 392.d0/429.d0
      return

  END FUNCTION bwint

!=======================================================================
   FUNCTION bwzeta(i1)
!=======================================================================
!
!  Computes the nuclear spin-orbit parameter and the
!  using the formula derived by blume and watson.
!
!
!----------------------------------------------------------------------
!
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: i1
     REAL(KIND=8) :: bwzeta

      PARAMETER(nt=100,kx=8,nwd=10)
      COMMON/blume/coefn2(4),coefnk(4),coefvk(4)
!
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
!
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      dimension ss(3)
!
!
      zeta = fine*z*quadr(i1,i1,-3)
      lb = l(i1)
      do i = 1,nwf
        if (i .eq. i1) cycle
        la = l(i)
        zeta = zeta -sum(i)*sn(i1, i, i1, i, 0)
        if (sum(i) .ne. 4*l(i)+2) go to 10
          call bwint(la,lb)
          ke1 = 2
          if (la .ne. lb) ke1 = iabs(la-lb)
          ip = 0
          do k = ke1,la+lb,2
            ip = ip+1
            zeta = zeta+coefn2(ip)*sn(i1, i, i, i1, k-2)
     :                 +coefnk(ip)*sn(i, i1, i1, i, k)
     :                 +coefvk(ip)*(vk(i1,i,i,i1,k-1)-vk(i,i1,i1,i,k-1))
          end do
      end do
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
      END SUBROUTINE bwzeta
!=======================================================================
        subroutine bxvpw(c,b,v,y,w)
!=======================================================================
!   Computes y = c* b * v + y   where b is a symmetric, banded matrix
!   and v, w are vectors
!
!   Written by C. F. Fischer
!------------------------------------------------------------------------
!
!
!   on entry
!   --------
!     nt      the leading dimension of arrays.
!     k       the number of diagonals
!     n       the order of the matrix
!     c,ic    coefficients
!     b       the symmetric, banded matrix in column storage mode
!     v       vector
!     w       working array
!
!   on exit
!   -------
!     y       y = c*B*v +y
!-----------------------------------------------------------------------
!
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: i,j,k
     REAL(KIND=8), INTENT(IN), DIMENSION(nt,:) :: b
     REAL(KIND=8), INTENT(IN), DIMENSION(:) :: v,y
     REAL(KIND=8), INTENT(INOUT), DIMENSION(:) :: y

!
! ...   initialize the w array
!
	 w = 0.d0
!
!      .. contribution from sub-diagonals
!
       do i=1,k
         do j=k-i+1,n
	   w(j) = w(j)  +b(j,i)*v(j-k+i)
         end do
       end do
!
!      .. contribution from super-diagonals
!
       do i=1,k-1
         do j=1,n-k+i
	   w(j) = w(j) + b(j+k-i,i)*v(j+k-i)
         end do
       end do
!
       if ( c .ne. 1.d0) then
           y = y + c*w
       else
           y = y + w
       end if
       END SUBROUTINE bxvpw
      
!=======================================================================
   FUNCTION ca(l,k) 
!=======================================================================
! Computes the direct average interaction
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l, k
    REAL(KIND=8) :: ca

      ca = rme(l,l,k)**2

    END FUNCTION ca

!=======================================================================
   FUNCTION cb(l,lp,k) 
!=======================================================================
!  Compute the average exchange interaction
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN): l,lp,k
    REAL(KIND=8):: cb

      cb = rme(l,lp,k)**2/(2*(2*l+1)*(2*lp+1))

    END FUNCTION cb
!=======================================================================
   SUBROUTINE data` 
!=======================================================================
!  Determines data about the problem to be solved:
!   . number of congifugration states (ncfg)
!   . number and types of orbitals
!   . data associated with the energy expression`
!----------------------------------------------------------------------
!
    IMPLICIT NONE

      PARAMETER(nt=100,kx=8,nwd=10)
*
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
      CHARACTER config*50,atom*6,term*6,ans*1,string*50
      COMMON /LABEL/config,el(nwd),atom,term
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      LOGICAL fail,omit,rel,all,trace,first,strong,done,orderd
      COMMON /TEST/fail,omit,rel,all,trace
*
      INTEGER ind(nwd)
      CHARACTER*3 el,el1,el2,elcsd(18)
      CHARACTER*1 aster,w
      data aster/'*'/
*
    7 format(a3,f6.0,i3,i3,f3.1)
*
*  *****  read 'atom' card
*
INPUT: do
     write(iscw,'(a/a)') ' Enter ATOM, TERM, Z', &
               ' Examples: O,3P,8. or Oxygen,AV,8.'
     read(5,'(a50)') string
    i = index(string,',')
    if ( i == 0) then
      write(iscw,*)' ATOM, TERM, and Z must be separated by commas '
      cycle INPUT
    else 
      atom = string(1:i-1)
      j = index(string(i+1:),',')
      if ( j == 0) then
        write(iscw,*)' ATOM, TERM, and Z must be separated by commas '
        cycle INPUT
      end if
      term = string(i+1:i+j-1)
      read(string(i+j+1:), '(f3.0)') z

      !  input COMMON closed shells
      write(iscw,*)
      write(iscw,'(a,a)')' List the closed shells in fields indicated',
     :        ' (blank line if none)'
      write(iscw,'(a)') ' ... ... ... ... ... ... ... ... etc.'
      read(5,'(18(1X,A3))') (elcsd(i),i=1,18)

      ! input the configurations
      write(iscw,'(/a,a/a)')' Enter electrons outside closed shells ',
     : '(blank line if none)',' example: 2s(1)2p(3)'
      read(5,'(a)')  string
      call reform(string, config)

      ! determine the number of closed shells
 
      i = 0
      ss = 0.d0
      lmax = 0
      if (elcsd(i+1) .ne. '   ') then
         i = i+1
         el(i) = elcsd(i)
         j = 3
         if (el(i)(1:1) .ne. ' ') j = 2
         l(i) = lval(el(i)(j:j))
         n(i) = ichar(el(i)(j-1:j-1)) - ichar('1') + 1
         ifull = 2*(2*l(i)+1)
         sum(i) = ifull
         s(i) = ss + ifull/2
         ss = ss + ifull
         meth(i) = 1
         acc(i) = 0.d0
         ind(i) = 0
         if (iuf .ne. 0)  ind(i) = -1
         if( i .lt. 18) go to 12
         stop ' too many closed shells: max = 18'
      end if
      nclosd = i
*
*  *****  determine the other electrons
*
      maxorb = nclosd
      string = config
      j = 2
      i = 0
 16   if (string(j:j+2) .ne. '   ' ) then
*
*  --------- an electron has been found; is it a new one?
*
         i = i+1
         if (i .gt. (5)) stop ' too many shells: max= (5)'
         el1 = string(j:j+2)
         k = nclosd + 1
 17      if (k .le. maxorb) then
            if ( el(k) .ne. el1 ) then
               k = k+1
               if (k .gt. (nwd)) then
                  write(iscw,'(a,i4)')' too many electrons: max =',nwd
                  go to 2
               else
                  go to 17
               end if
            end if
         else
*
*  ------------  a new electron has been found; add it to the list
*
            maxorb = k
            el(maxorb) = el1
            read(string(j+4:j+7),'(f4.0)') sum(k)
         end if
         j = j+10
         if (j .lt. (50)) go to 16
      end if
*
*  -----  the list of electrons has been determined
*
      write(iscw,19) maxorb,(el(j),j=1,maxorb)
   19 format(/' there are ',i3,' orbitals as follows:'/(1x,18(1x,a3)))
      nwf = maxorb
      if (nit .lt. 0) nit=nwf
21    write(iscw,'(/a,a)') ' orbitals to be varied: ',
     :  'all/none/=i (last i)/comma delimited list/h'
      read(5, '(a)') string
      if (string(1:1) .eq. 'h' .or. string(1:1) .eq. 'h') then
       call help(1)
       go to 21
      else if (string(1:3) .eq. 'all' .or. string(1:3) .eq. 'all') then
         nit = nwf
      else if (string(1:4).eq.'none' .or. string(1:4).eq.'none') then
         nit = 0
      else if (index(string,'=') .ne. 0) then
         j = index(string,'=')
	 jj = index(string,' ')
	 if ( jj .eq. j+2) then
	    read(string(j+1:j+1),'(i1)') nit
	 else if (jj .eq. j+3) then
	    read(string(j+1:j+2),'(i2)') nit
	 else
	    write(iscw,'(a,a)') ' nit must be specified by one or two',
     :        'digits immediately following = sign (no blanks): Renter'
	    go to 21
	 end if
      else
         nit = 0
         j = 1
22       next = index(string(j:),',')
*
*        ***  search for last electron LABEL which need not be followed
*             by a comma
*
         if (next .eq. 0 .and. string(j:j+2) .ne. '   ')
     :       next = index(string(j+1:),' ') + 1
         if (next .ge. 1) then
            if (next .eq. 4) then
               el1 = string(j:j+2)
            else if (next .eq. 3) then
               el1 = ' '//string(j:j+1)
            else
               write(iscw,*)
     :              ' electron LABELs must be separated by commas;'
               write(iscw,*)' each LABEL must contain 2 or 3 CHARACTERs'
               go to 21
            end if
            call reord(el,el1,nwf,ierr)
            if (ierr .eq. 0) then
               nit = nit + 1
               j = j + next
               if (j .lt. 72) go to 22
            else
               write(iscw,*) ' case must match as well as position of',
     :                  ' imbedded blanks'
               write(iscw,*) ' for 3rd CHARACTER of LABEL to be blank',
     :                 ' follow blank with comma'
     
               go to 21
            end if
         end if
      end if
*
      ib = nwf - nit + 1
      if (nit .ne. 0) then
23       write(iscw,'(/a)') ' default electron parameters ? (y/n/h) '
         read(5,'(a)') ans
         if ( ans .eq. 'h' .or. ans .eq. 'h' ) then
            call help(2)
            go to 23
         end if
      else
         ans = 'y'
*
      end if
      if ( ans .ne. 'y' .and. ans .ne. 'y') write(iscw,'(a,a)')
     :   ' s, ind, meth, acc for non-closed shell electrons: '
      do 20 i = nclosd+1,nwf
         if ( ans .eq. 'y' .or. ans .eq. 'y') then
            s(i) = ss + (sum(i)-1.d0)/2.d0
            ss = ss + sum(i)
            meth(i) = 1
            acc(i) = 0.d0
            ind(i) = 0
            if (iuf .ne. 0)  ind(i) = -1
         else
            write(iscw,'(a,a)') el(i),':  '
            read(5,*) s(i),ind(i),meth(i),acc(i)
         end if
         j = 2
         if (el(i)(1:1) .eq. ' ') j = 3
         l(i) = lval(el(i)(j:j))
         n(i) = ichar(el(i)(j-1:j-1)) - ichar('1') + 1
 20   continue
      if (nwf .eq. 1) meth(1) = 2
*
*  *****  define all orbitals in the configuration to be orthogonal
*
       do 34 i = 1,nwf
         e(i,i) = 0.d0
	 lmax = max0(lmax,l(i))
         do 33 j = 1,i-1
            e(i,j) = 0.d0
            if (l(i) .eq. l(j)) e(i,j) = 1.d-5
            e(j,i) = e(i,j)
   33       continue
   34 continue
      kmax = 2*lmax
      ib = nwf - nit + 1
      no = ns
      strong = .false.
      write (3,62) atom,term,z,(el(i),int(sum(i)),i=1,nclosd)
62    format(1h1///9x,33hhartree-fock wave functions for  ,2a6,4h z =,
     : f5.1//14x,'core =',5(1x,a3,'(',i4,')')/(20x,5(1x,a3,'(',i4,')')))
      write (3,'(5x,a15,a50)')  'configuration =',config
68    continue
      write (3,71)
71    format(//9x,10hinput data/9x,10h----- ----//13x,13hwave function,
     :   11h  procedure/17x,22hnl  sigma meth acc opt///)
      do 79 i = 1,nwf
         write (3,78) i,el(i),n(i),l(i),s(i),meth(i),acc(i),ind(i)
78       format(i8, 2x,a3,2i3,f7.1,i4,f4.1,i4)
79    continue
      omit = .not. strong
*
      call array
      call enexpr(term, done)
      if (.not. done) then
*
*  ---  case needs additional data
*
         write(iscw,85)
85       format(/' The program could not derive the energy expression'/
     :           ' Select one of the following options and enter:'/
     :           '    1  re-enter the term and configuration'/
     :           '    2  enter the deviations from eav as input'/
     :           '    3  stop'/)
         read(5,*) iselec
         go to (5,86,99) iselec
86       write(iscw,83)
83       format(/' Input data for deviations from the average energy'/
     :   ' first fk integrals, then gk integrals in indicated format'/
     :   '  cc.ccccccccccfkk(el1,el2)  - terminate each list with an *',
     :   ' in the f column')
         first = .true.
*
*  *****  read 'fk' and 'gk' cards, omitting the header if a file
*
82       read (5,84) cfg,w,kfg,el1,el2
84       format(f14.8,a1,i2,1x,a3,1x,a3)
         if ( w .ne. aster ) then
            call eptr(el,el1, ifg, *99)
            call eptr(el,el2, jfg, *99)
            call add(cfg,kfg,ifg,jfg,first)
            call add(cfg,kfg,jfg,ifg,first)
            go to 82
         else if (first) then
            first = .false.
            go to 82
         end if
      end if
*
*  *****  compute the initial array and initial radial functions
*
      call initsp
      call wavefn(ind)
*
*      ... define an order for the functions to be iterated
*
      do 90 jp = 1,nwf
         iord(jp) = jp
90    continue
91    orderd = .true.
      do 92 jp = ib,nwf-1
         n1 = n(iord(jp))
         l1 = l(iord(jp))
         n2 = n(iord(jp+1))
         l2 = l(iord(jp+1))
         if (n1.gt.n2 .or. (n1.eq.n2 .and. l1.gt.l2)) then
            itemp = iord(jp)
            iord(jp) = iord(jp+1)
            iord(jp+1) = itemp
            orderd = .false.
         end if
92    continue
      if (.not. orderd) go to 91
      return
99    stop
      end
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
!=======================================================================
   SUBROUTINE dev(iel, l, q, i, done) 
!=======================================================================
!
!     add the deviations to the average energy for a partially filled
!       p- or d- shell
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN):: iel, l, q, i
    LOGICAL, INTENT(OUT) :: done 
      INTEGER f2pp(6), f2dd(45), f4dd(45)
      data    f2pp/-3,3,12,-9,0,6/
      data    f2dd/-58,77,50,-13,140,
!             ... d3 coefficients
     :       -93,42,-12,-57,123,105,69,-12,
!             ... d4 coefficients
     :        -105,-69,-24,66,12,39,21,57,
     :        -51,30,48,84,219,111,210,138,
!             ... d5 coefficients
     :        -175,-85,23,-22,-112,-76,-58,167,
     :        23,-85,59,140,104,86,320,113/
      data    f4dd/5,-70,15,50,140,
!             ... d3 coefficients
     :        -30,-105,30,55,-45,105,-15,30,
    
!             ... d4 coefficients
     :        -105,15,-10,45,-30,-45,70,-55,
     :        75,135,20, 0,30,-15,210,-30,
!             ... d5 coefficients
     :        -175,-50,-40,-85,35,50,110,-15,
     :        -5,125,-25,140,20,-40,-100,-55/
     
      done = .true.
      n = q
      if (n .gt. 2*l+1) n = 4*l+2 - n
      if (n .gt. 1) then
         if (l .eq. 1) then
            call add(2*f2pp(i)/25.d0,2,iel,iel,.true.)
         else if (l .eq. 2) then
            i = i-6
            call add(2*f2dd(i)/441.d0,2,iel,iel,.true.)
            call add(2*f4dd(i)/441.d0,4,iel,iel,.true.)
	 else
	    done = .false.
         end if
      end if
      return
      END SUBROUTINE dev
!=======================================================================
   FUNCTION ekin(i,ii,rel) 
!=======================================================================
!
!       returns the value of the integral of
!
!         (2/r)p (y p  + x )
!               j  i i    i
!
!   integrated with respect to r.
!
!----------------------------------------------------------------------
!
    IMPLICIT NONE

!
      double precision function ekin(i,ii,rel)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(nt=100,kx=8,nwd=10)
!
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
!
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
!=======================================================================
   SUBROUTINE energy(etotal) 
!=======================================================================
!
!       determines the total energy
!
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    REAL(KIND=8), INTENT(OUT) :: energy

      PARAMETER(nt=100,kx=8,nwd=10)
      CHARACTER config*50,el*3,atom*6,term*6
      COMMON /LABEL/config,el(nwd),atom,term
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      LOGICAL fail,omit,rel,all,trace
      COMMON /TEST/fail,omit,rel,all,trace
      DOUBLE PRECISION ek(nwd)
*
*  *****  compute kinetic energy if necessary
*
      do i = 1,nwf
        ek(i) = -0.5d0*hl(el,i,i,rel)
      end do
*
      etotal = 0.d0
      do i = 1,nwf
         etotal = etotal + sum(i)*(ek(i))
	 print *, 'etotal', ek(i), etotal
         do j = 1,i
            do k = 0,2*min0(l(i),l(j)),2
               c = a(i,j,k)*sum(i)
               if (i .eq. j) c = c/2.d0
               if (abs(c).ne.0.d0) etotal = etotal + c*fky(i,j,k,rel)
	       print *, i,j,c,fky(i,j,k,rel),etotal
            end do
         end do
         do j = 1,i-1
           do k = abs(l(i)-l(j)),l(i)+l(j),2
             c = b(i,j,k)*sum(i)
             if (abs(c).ne.0.d0) etotal=etotal+c*gky(i,j,k,rel)
	       print *, i,j,c,gky(i,j,k,rel),etotal
           end do
         end do
       end do
     END SUBROUTINE energy
*
*     ------------------------------------------------------------------
*                       e n e x p r
*     ------------------------------------------------------------------
*
*     determine the deviations to the average energy for the following:
*        i) an open p- or d-shell
*       ii) a single electron or hole, any l
*      iii) an s-electron and a single electron, any l
*       iv) an s-electron and an open p- or d-shell
*        v) an open p-shell and a single electron, any l
*
      SUBROUTINE enexpr(term, done)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(nt=100,kx=8,nwd=10)
      CHARACTER term*6,sl*2,senor*1,psl*2,slm*2,slp*2
      LOGICAL done
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      INTEGER sumtab(5),partab(11),ptrtab(11),ltab(54),
     :   nos(2),plval(11),pacval,sp,ps1,ps2
*
*     ... fint, gint1, and gint2 are coefficients of polynomials
*         in l, tabulated by slater,
*
*
      INTEGER fint(3,54),gint1(3,54),gint2(3,54)
      CHARACTER*1 parch(11)
*
*     ... coefficients of f2 integrals for p(n)l(1) configurations
*
      data fint/2,-1,0,-4,-4,3,2,5,3,2,-1,0,-4,-4,3,2,5,3,
     :   -2,1,0,4,4,-3,-2,-5,-3,-2,1,0,4,4,-3,-2,-5,-3,
     :   4,-2,0,-2,-11,6,-4,-4,15,-2,7,15,4,10,6,0,0,0,
     :   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     :   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     :   2,-1,0,-4,-4,3,2,5,3,2,-1,0,-4,-4,3,2,5,3,
     :   -4,2,0,2,11,-6,4,4,-15,2,-7,-15,-4,-10,-6,0,0,0,
     :   -2,1,0,4,4,-3,-2,-5,-3,-2,1,0,4,4,-3,-2,-5,-3/
*
*     ... coefficients of g(l-1) integrals
*
      data gint1/-10,5,0,2,11,-6,2,-1,-6,14,-7,0,2,-13,6,2,-1,6,
     :   -8,4,0,-8,4,0,4,10,0,10,-5,0,10,-5,0,4,-8,0,
     :   -8,4,0,-2,13,-6,2,11,-12,4,4,-12,4,-2,-12,0,0,0,-6,3,0,10,-5,0,
     :   -6,3,0,-6,3,0,-4,8,3,0,12,3,6,9,0,-6,3,0,6,21,12,12,12,-27,
     :   12,-6,27,6,-15,-12,-6,3,0,0,6,-3,0,0,-3,6,-3,0,0,-6,3,12,6,3,
     :   -4,2,0,-4,2,0,-4,2,0,-4,2,0,14,11,-9,14,-7,-9,-4,-2,0,-4,2,0,
     :   -2,7,3,2,11,3,8,8,0,0,0,0,-2,1,0,-2,1,0,-2,1,0,-2,1,0,-2,1,0,
     :   22,13,0/
*
*     ... coefficients of g(l+1) integrals
*
      data gint2/2,5,-3,2,-7,-15,-10,-25,-15,2,5,9,2,17,21,14,35,21,
     :   4,-2,-6,-8,-20,-12,-8,-20,-12,4,16,12,10,25,15,10,25,15,
     :   4,10,0,4,4,-12,2,-7,-21,-2,-17,-21,-8,-20,-12,0,0,0,
     :   -6,-15,-9,10,25,15,6,3,-3,0,-12,-9,-4,-16,-9,-6,-15,-9,
     :   -6,-15,-9,6,27,9,12,30,-9,12,12,-27,6,-9,-27,-6,-15,-9,
     :   0,0,-3,0,-6,-9,-6,-15,-9,12,18,9,0,6,9,6,15,9,
     :   -4,-10,-6,-4,-10,-6,-4,-10,-6,14,35,12,14,17,-6,-4,-10,-6,
     :   8,8,0,2,-7,-6,-2,-11,-6,-4,-10,-6,-4,-10,-6,0,0,0,
     :   -2,-5,-3,-2,-5,-3,-2,-5,-3,22,31,9,-2,-5,-3,-2,-5,-3/
*
*     ... encoded term value -- s = ltab/10
*                         lterm = l + (ltab mod 10 - 5)
*         example:  ltab = 36 with l = 2  is 3f
*
      data ltab/36,35,34,16,15,14,46,45,44,26,25,24,27,26,25,24,
     :   23,25,55,35,37,36,35,34,33,17,16,15,14,13,36,35,34,16,15,14,
     :   46,45,44,26,25,24,27,26,25,24,23,25,36,35,34,16,15,14/
      data sumtab/1,4,7,10,11/
      data partab/2,3,1,1,4,2,2,3,1,1,2/
      data ptrtab/6,12,17,18,20,30,36,42,47,48,54/
      data plval/1,1,2,0,0,2,1,1,2,0,1/
      data parch/'p','p','d','s','s','d','p','p','d','s','p'/
*
      ip = 1
 1    if (term(ip:ip) .eq. ' ') then
         ip = ip+1
         go to 1
      end if
      sl = term(ip:ip+1)
      senor = ' '
      if (ip.le.4) senor = term(ip+2:ip+2)
*
*   ---  convert lowercase l symbol to uppercase
*
      if (sl(2:2).gt.'a' .and. sl(2:2).lt.'z')
     :    sl(2:2) = char(ichar(sl(2:2)) + ichar('a') - ichar('a'))
*
*  ---  determine if fk or gk data needs to be input
*
      il = 0
      is = 0
      j = 1
      do 2 i = nclosd+1, nwf
         if (sum(i) .ne. 4*l(i)+2 .and. sum(i) .ne. 0.d0) then
            if (j.gt.2) then
               if (sl .ne. 'av' .and. sl .ne. 'av') then
                  done=.false.
               else
                  done = .true.
               end if
               return
            endif
            nos(j) = i
            j = j + 1
            if (l(i) .eq. 0 .and. is .eq. 0) then
               is = is + 1
               iis = i
            else
               il = il + 1
               iil = i
            end if
         end if
 2    continue
      if (sl .ne. 'av' .and. sl .ne. 'av' .and. is+il.ne.0) then
         done = .false.
         c = 0.d0
         if (is+il .le. 2 .and. il .le. 1) then
            if (is .eq. 0 .and. il .eq. 1) then
3              call looktm(l(iil),sl,senor,sum(iil),ip,nsl)
	       if (nsl .gt. 1) then
		  write(iscw,*)' ambiguous term: enter seniority'
		  read (5,'(a1)') senor
		  go to 3
	       end if
               call dev(iil,l(iil),sum(iil),ip,done)
            else if (is .eq. 1 .and. il .eq. 1) then
               slm = sl
	       slp = sl
               slm(1:1) = char(ichar(slm(1:1)) - 1)
               slp(1:1) = char(ichar(slp(1:1)) + 1)
	       call looktm(l(iil),slm,senor,sum(iil),ipm,nslm)
	       call looktm(l(iil),slp,senor,sum(iil),ipp,nslp)
	       if (nslm+nslp .eq. 0) then
	          done = .false.
	          return
	       else if (nslm .eq. 1 .and. nslp .eq. 0) then
		  sl = slm
		  ip = ipm
	       else if (nslm .eq. 0 .and. nslp .eq. 1) then
		  sl = slp
		  ip = ipp
	       else if (nslm .eq. 1 .and. nslp .eq. 1) then
4                 write(iscw,'(a,a3,a,a3)')
     :                 ' ambiguous l**n term: enter',slm, ' or ',slp 
	          read(5,'(a2)') sl
		  if (sl .eq. slm) then
		     ip = ipm
		  else if (sl .eq. slp) then
		     ip = ipp
		  else
	             write(iscw,*) ' term not allowed: re-enter'
		     go to 4
		  end if
	       else 
5                 write(iscw,'(a,a)') ' ambiguous l**n parent term:',
     :                   'enter term and seniority'
	          read(5,'(a2,a1)') sl, senor
	          if (senor .eq. ' ') then
			write(iscw,*) 'seniority is needed'
			go to 5
	          end if
	          call looktm(l(iil),sl,senor,sum(iil),ip,nsl)
	          if (nsl .ne. 1) then
		     write(iscw,'(a,a3,a,a3,a)') ' allowed terms are ',
     :                    slm, ' or ', slp,' plus seniority'
		     go to 5
	          end if
	       end if
               call dev(iil,l(iil),sum(iil),ip,done)
               if (done ) then
                  sp = ichar(sl(1:1)) - ichar('0')
                  csp = (sp - 1)/2.
		  if (sl .eq. slm) then
                     c = -csp/(2*l(iil)+1)
                  else
                     c = (csp + 1)/(2*l(iil)+1)
                  end if
                  call add(c,l(iil),iis,iil,.false.)
                  call add(c,l(iil),iil,iis,.false.)
               end if
            else if (is .eq. 1 .and. il .eq. 0) then
               done = .true.
            end if
         else
            if (((l(nos(1)).eq.1).and.(sum(nos(2)).eq.1.d0)).or.
     :          ((l(nos(2)).eq.1).and.(sum(nos(1)).eq.1.d0))) then
               if (l(nos(1)).eq.1.and.sum(nos(2)).eq.1.d0) then
                  isump=sum(nos(1))
                  np = nos(1)
                  nl = nos(2)
               else
                  isump = sum(nos(2))
                  np = nos(2)
                  nl = nos(1)
               endif
               sp=ichar(sl(1:1))-ichar('0')
               lp=lval(sl(2:2))
               ps1=sp+1
               ps2=sp-1
               if (isump.eq.1) then
                  iptr1=1
               else
                  iptr1=sumtab(isump-1)+1
               end if
               iptr2=sumtab(isump)
               nomach=0
               call lookup(partab,iptr1,iptr2,ind,nomach,ps1)
               call lookup(partab,iptr1,iptr2,ind,nomach,ps2)
               psl(1:1)=char(partab(ind)+ ichar('0'))
               psl(2:2)=parch(ind)
               if (nomach.gt.1) then
                  write(iscw,*)' ambiguous parent case'
 10               write(iscw,*)' enter the sl term for p(n) subshell'
                  read (5,'(a)')psl
                  if (psl(2:2).gt.'a'.and.psl(2:2).lt.'z')
     :                   psl(2:2)=char(ichar(psl(2:2))+ichar('a')
     :                                                -ichar('a'))
                  ps1=ichar(psl(1:1))-ichar('0')
                  ps2=lval(psl(2:2))
                  call lookup(plval,iptr1,iptr2,ind,nomach,ps2)
                  if ((nomach.ne.1).and.(partab(ind).ne.ps1))
     :                    go to 10
               end if
               if (isump.eq.1) then
                  iptr1=1
               else
                  iptr1=ptrtab(ind-1)+1
               end if
               iptr2 = ptrtab(ind)
               lv=l(nl)
               pacval=sp*10+lp-lv+5
               nomach=0
               call lookup(ltab,iptr1,iptr2,ind,nomach,pacval)
               if (nomach.ne.1) then
                  done=.false.
                  return
               endif
               val1=((fint(1,ind)*lv+fint(2,ind))*lv+fint(3,ind))
     :                  /(5.d0*(2*lv-1)*(2*lv+3))
               val2=((gint1(1,ind)*lv+gint1(2,ind))*lv+gint1(3,ind))
     :                  /(2.d0*(2*lv+1)*(2*lv-1)**2)
               val3=((gint2(1,ind)*lv+gint2(2,ind))*lv+gint2(3,ind))
     :                  /(2.d0*(2*lv+1)*(2*lv+3)**2)
*
*     ...  add contributions from between p-subshell and l-electron
*
               call add(val1,2,np,nl,.true.)
               call add(val1,2,nl,np,.true.)
               call add(val2,lv-1,np,nl,.false.)
               call add(val2,lv-1,nl,np,.false.)
               call add(val3,lv+1,np,nl,.false.)
               call add(val3,lv+1,nl,np,.false.)
*
*     ... add deviations for p-subshell
*
               call looktm(1,psl,' ',sum(np),ip,nsl)
               call dev(np,1,sum(np),ip,done)
            else
               done = .false.
            end if
         end if
      else
         done = .true.
      end if
      return
      end
*
*     ------------------------------------------------------------------
*                       e p t r
*     ------------------------------------------------------------------
*
*       determines the position of the electron in the electron list
*
      SUBROUTINE eptr(el,elsymb, iel, *)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER (nt=100,kx=8,nwd=10)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      CHARACTER el(*)*3, elsymb*3, bl*3
      data bl/'   '/
*
* ***** search electron list for lsymb
*
      if ( elsymb .eq. bl ) then
         iel = 0
         return
      endif
      do 10 i=1,nwf
         if (el(i) .eq. elsymb ) then
            iel = i
            return
         endif
10    continue
      iel = -1
      write(iscw,20) elsymb
20    format(/10x,a3,' not found in electron list')
      return 1
      end
*
*     -----------------------------------------------------------------
*           f a c t r l
*     -----------------------------------------------------------------
*
*
      SUBROUTINE factrl(nfact)
*
*      gam(i) = log( gamma(i-1) ), where gamma(i) = factorial i-1
*
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
*
      COMMON/fact/gam(100)
      data zero,one,two/0.d0,1.d0,2.d0/
*
      gamma=one
      gam(1) = zero
      do 1 i=1,nfact-1
         gamma=i*gamma
         gam(i+1) = dlog(gamma)
    1 continue
      do 20 i = nfact+1,(100)
         x = i-1
         gam(i) = gam(i-1) + dlog(x)
   20 continue
      return
      end

      subroutine ginvit(a,lda,n,c,ipvt,x,ierr)
      implicit real*8(a-h,o-z)
      dimension ipvt(lda),a(lda,*),c(*),x(*)
      data eps/1d-14/
*
*     Generalized inverse Iteration problem
*               Ac = 0
*     where
*             A = A* + k^2*B
*             A is dense matrix
*
*    lda  -  leading dimension of a
*    n    -  size of the matrix
*    c    -  initial guess of the eigenvector on input,
*            final solution on output
*    ipvt -  space for pivot information
*    x    -  working vector
*    ierr -  0 if converged, number of iterations if 30.
*
*    call for LU decomposition of A
*     print *, 'Entering ginvit with n,c =',n,(c(i),i=1,n)
      call dgefa(a,lda,n,ipvt,info)
CTSS  call sgefa(a,lda,n,ipvt,info)
CUNIC call sgefa(a,lda,n,ipvt,info)
      if(info.ne.0)then
	 write(0,*) 'Call to dgefa produced zero in position ',info
CTSS	 write(0,*) 'Call to sgefa produced zero in position ',info
CUNIC    write(0,*) 'Call to sgefa produced zero in position ',info
	 a(info, info) = 1.d-14
      endif

c   Inverse Iteration starts from here

      ierr = 0
      it = 0
      aold = 0
 200  it = it + 1
      call dgesl(a,lda,n,ipvt,c,0)
CTSS  call dgesl(a,lda,n,ipvt,c,0)
CUNIC call dgesl(a,lda,n,ipvt,c,0)
      anorm = 0.
      amaxx = 0.
      do 30 j =1,n
	 x(j) = c(j)
	 absx = dabs(x(j))
	 anorm = anorm + absx
	 if(absx.gt.amaxx) then
	    amaxx = absx
	    jmax = j
	 end if
 30   continue
C     write(6,*) it,' Maximum is x(',jmax,') = ',x(jmax)
      del = anorm - aold
      aold = anorm
      err = dabs(del/anorm)
      scale = 1.d0/x(jmax)
      do 40 i = 1,n
	 x(i) = x(i)*scale
	 c(i) = x(i)
 40   continue
      if(err.gt.eps .and. it .lt. 30) then
c        do 15 i = 1,n
c            x(i) = 0d0
c            do 17 j = 1,n
c               x(i) = x(i) + b(i,j)*c(j)
c17          continue
c15      continue
c        do 18 i = 1,n
c           c(i) = x(i)
c18      continue
         go to 200
      else
         do 50 i  = 1,n
           c(i)   = x(i)
 50      continue
	 if (it .eq. 30) then
	    ierr = 30
	 end if
      endif
      return
      end
*     ==================================================================
*              G R A D
*     ==================================================================
*
*  *****  THE GRAD FUNCTION SUBPROGRAM COMPUTES THE FOLLOWING DIRECTLY
*  *****         <P(J)^D + L(I)/R ^P(I)> WITH L(I) > L(J)
*
      DOUBLE PRECISION  FUNCTION GRAD(I,J,ier)
      IMPLICIT DOUBLE PRECISION (a-h,o,p,r-z)
      PARAMETER (nt=100,kx=8,nwd=10)
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rk(nt,nt,kx,kx),lh,krk
*
      LOGICAL ixch
      DIMENSION g(nt,kx),a(nt,kx)
      
      IF ( IABS(L(I) - L(J)) .NE. 1) THEN
	ier = 1
	grad = 0.0d0
      else
      LL = MAX0(L(I),L(J))
      i1 = i 
      i2 = j
      ixch = .false.
      IF ( L(I) .lt. L(J) ) THEN
	i1 = j
	i2 = i
	ixch = .true.
      end if
*
* ...  form the sums
*
      last = min0(max(i1),max(i2))
      grad = 0.d0
      do 100 jj = 2,last
         grad = grad + p(jj,i1)*p(jj,i2)*rm1(jj,ks)
  100 continue
      grad = ll*grad
      do 110 jp = 2,ks
*         .. j' < jj
        do 110 jj = 2,last-jp+1
          grad = grad +
     :     g(jj,ks-jp+1)*(p(jj,i2)*p(jj-jp+1,i1)-
     :                    p(jj,i1)*p(jj-jp+1,i2)) +
     :     ll*rm1(jj,ks-jp+1)*(p(jj,i2)*p(jj-jp+1,i1)+
     :                         p(jj,i1)*p(jj-jp+1,i2)) 
  110     continue

      IF ( ixch ) GRAD = - GRAD
      end if
      END
*               g r a n g e
*     ------------------------------------------------------------------
*
*       controls the calculation of off-diagonal energy parameters.
*   it searches for all pairs (i,j) which are constrained through  an
*   orthogonality requirement.   eq. (7-10) is used to calculate the
*   PARAMETER.
*
      SUBROUTINE grange
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(nt=100,kx=8,nwd=10)
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
      LOGICAL fail,omit,rel,all,trace
      COMMON /TEST/fail,omit,rel,all,trace
*
*  *****  rotate pairs connected by orthogonality but not when one of
*         the orbitals is simultaneously orthogonal to a non-orthogonal
*         pair
*
      do 1 i = ib,nwf-1
         do 2 j = i+1,nwf
            if (dabs(e(i,j)) .gt. 1.d-10) call rotate(i,j)
2        continue
1     continue
      print *, 'sumi, sumj', sum(1), sum(2)
*
*  *****   compute off-diagonal energy parameters
*
      do 10 i = max0(2,ib),nwf
         do 12 j = 1,i-1
           if (dabs(e(i,j)) .gt. 1.d-10) then
             if ( j .lt. ib) then
               e(i,j) = hl(el,i,j,rel) - ekin(i,j,rel)
               e(j,i) = 0.d0
             else if (sum(i) .eq. sum(j)) then
               c=hl(el,i,j,rel)-(ekin(i,j,rel)+ekin(j,i,rel))/2.d0
               e(i,j) = c
               e(j,i) = c
             else
               res = 0.d0
               do 14 ii = 1,nwf
                 if (ii.eq.i .or. ii.eq.j) then
                   do 22 k = 0,2*l(i),2
                     if(ii.eq.i) then
                       c= a(i,i,k)-a(j,i,k)-b(j,i,k)
                       if (dabs(c).gt.1.d-10)
     :                          res = res + c*rky(i,i,i,j,k,rel)
		       print *, i,j,k,c,res
                     else if (ii.eq.j) then
                       c= a(j,j,k)-a(i,j,k)-b(i,j,k)
                       if (dabs(c).gt.1.d-10)
     :                         res = res - c*rky(j,j,j,i,k,rel)
		       print *, i,j,k,c,res
                     end if
22                continue
                else
                  do 24 k = 0,2*min0(l(i),l(ii)),2
                    c = a(i,ii,k) - a(j,ii,k)
                    if (dabs(c).gt.1.d-10)
     :                          res = res + c*rky(i,ii,j,ii,k,rel)
                    kk = abs(l(i)-l(ii)) + k
                    c = b(i,ii,kk) - b(j,ii,kk)
                    if(dabs(c) .gt. 1.d-10)
     :                            res = res + c*rky(i,ii,ii,j,kk,rel)
		       print *, i,j,k,c,res
24                continue
                end if
14            continue
              e(i,j) = 2.d0*sum(j)*res/(sum(i)-sum(j))
              e(j,i) = sum(i)*e(i,j)/sum(j)
            end if
          end if
          if (dabs(e(i,j)) .gt. 1.d-10) then
            write(6,35) el(i),el(j),e(i,j),el(j),el(i),e(j,i)
35          format(7x,2(3x,'e(',2a3,') =',f12.5))
          end if
 12     continue
 10   continue
      return
      end
     
*
*     ------------------------------------------------------------------
*               h e l p
*     ------------------------------------------------------------------
*
*       provide help information about the data requested
*
      SUBROUTINE help(case)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
      INTEGER case
*
      go to (10,20,30,40,50,60,70,80,90,100) case
*
10    write(iscw,11)
*  ***** which orbitals varied?
11    format(//1x,'response  all will vary all orbitals'/
     :         11x,'none will not vary any orbitals'/
     :         11x,'=n (INTEGER n), will vary last n'/
     :         11x,'with comma delimited list will vary',
     :         1x,'only the orbitals in the list'//)
      return
     
20    write(iscw,21)
*  ***** default electron parameters ?
21    format(//1x,
     :  'response  n  will prompt the user for:'//
     :   7x,'s     : screening parameter  (real number) '/7x,
     :  'ind   : indicator specifying the type of initial estimate'/
     :  15x,'0 - screened hydrogenic functions'/
     :  14x,'-1 - search for functions in wavefunction file; if'/
     :  20x,'not present use screened hydrogenic.'/
     :  7x,'meth  : method for solving the differential equation'/
     :  15x,'1 - method 1 solves the boundary value problem for an'
     :  /20x,'acceptable solution which need not be normalized.'/
     :  15x,'2 - method 2 solves the boundary value problem for '/
     :  20x,'an acceptable solution which is normalized to first'/
     :  20x,'order. if the exchange function is identically zero the '/
     :  20x,'program will automatically select method 2.'/
     :  15x,'3 - method 3 is similar to method 1 but omits all checks'/
     :  20x,'for acceptability.'/
     :  7x,'acc   : inital accelerating factor '/
     :  18x,'( real number such that 0 .le. acc .lt. 1 )'//)
      return
     
30    write(iscw,31)
*  ***** default values (no,strong) ?
31    format(//1x,
     :  'response  y will set default values as follows'//
     :  10x,'no - maximum number of points in the range of the '/
     :  15x,'function is set to 200.'/
     :  10x,'strong - is LOGICALly set to .false.'//
     :  1x,'response n  will prompt the user for'//
     :  10x,'no - maximum number of points in the range of '/
     :  15x,'the function which should be a positive INTEGER'/
     :  15x,'from 160 for a small atom to 220 for a large atom.'/
     :  10x,'strong - may be set to .true. or .false by user.'//)
      return
     
40    write(iscw,41)
*  ***** default values for remaining parameters ?
41    format(//1x,
     :  1x,'response  y --sets the following default values'//
     :  10x,'print=.false.'/10x,'scftol=1.d-8'/10x,'nscf=12'/
     :  10x,'ic=2 + (nwf + 1 - ib)/4'/10x,'trace=.false.'//
     :  1x,'response  n --prompts user for new parameter values.'//)
      return
     
50    write(iscw,51)
*  ***** default values for print, scftol ?
51    format(//1x,
     :  1x,'response  y --default value for print is .false. thus'/
     :  16x,'radial functions are not printed.'//
     :  16x,'scftol -the initial value of the parameter defining the'/
     :  16x,'self-consistency tolerance for radial functions is set'/
     :  16x,'to a default value of  1.d-8 .'//
     :  1x,'response  n --prompts user for new values of print and ',
     :  'scftol .'//)
      return
     
60    write(iscw,61)
*  ***** default values for nscf,  ic ?
61    format(//1x,
     :  1x,'response  y --nscf, the maximum number of cycles for'/
     :  16x,'the scf process is set to default value of 12 and'/
     :  'ic is set to 2 + (nwf + 1 -ib)/4 .'//
     :  1x,'response  n --user prompted for new nscf and ic values'//)
      
     
70    write(iscw,71)
*  ***** default values for trace ?
71    format(//1x,
     :  'response  y --sets trace to default value of .false. thus a'/
     :  16x,'trace of energy adjustment will  not  be printed.'//
     :  1x,'response  n --prompts the user for the new value of trace.'/
     :  16x,'if .true. a trace will be printed showing the energy ',
     :  'adjustment'/16x,'process used by methd1 for finding an',
     :  ' acceptable solution'/16x,'with the correct number of',
     :  ' nodes.'//)
      return
     
80    write(iscw,81)
*  ***** additional parameters ?
81    format(//1x,
     :  'response  y --additional values may be computed for '/
     :  16x,'slater or magnetic integrals'/
     :  16x,'expectation values of r**k'/
     :  /16x,'electron density at the nucleus'/16x,'spin-orbit ',
     :  'parameter'/16x,'transition integrals'//
     :  1x,'response  n --additional parameter computation is ',
     :  'skipped'//)
      return
     
90    write(iscw,91)
*  ***** do you wish to continue along the sequence ?
91    format(//1x,
     :  'response  y --sequence is continued'//
     :  1x,'response  n --current case is ended, but new case may be'/
     :  16x,'started'//)
      return
     
100   write(iscw,101)
*  ***** do you wish to continue ?
101   format(//1x,
     :  'response  y --prompts user for additional iterations and'/
     :  16x,'new ic.  then performs additional nscf self-consistent'/
     :  16x,'field iterations.'//
     :  1x,'response  n --terminates the calculation.'//)
      return
      end
*
*     ------------------------------------------------------------------
*              H L
*     ------------------------------------------------------------------
*
*       Returns the value of <i|L|j> .. assumes symmetry 
*
      DOUBLE PRECISION  FUNCTION HL(EL,i,j,REL)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER (nt=100,kx=8,nwd=10)
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rk(nt,nt,kx,kx),lh,krk
*
      PARAMETER(IWRITE=6)
      CHARACTER EL(*)*3
      LOGICAL REL
*
      IF (IABS(L(i)-L(j)) .NE. 0) GO TO 999 
      if (lh .ne. l(i)) then
         call hlm(z,l(i),nt,kx,ks,ns,db2,rm1,rm2,h0)
	 lh = l(i)
      end if
      hl = 0.d0
      last = min0(max(i),max(j))
      do 100 k = 2,last
	hl = hl + p(k,i)*p(k,j)*h0(k,ks)
  100 continue
      if (i .eq. j) then
        do 110 kp = 2,ks
*         .. k'<k
	do   110 k = kp+1,last
	     hl = hl + 2.d0*p(k,i)*p(k-kp+1,j)*h0(k,ks-kp+1)
  110   continue
       else
        do 120 kp = 2,ks
*         .. k'<k
	  do 120 k = kp+1,last
	     hl = hl + 
     :         (p(k,i)*p(k-kp+1,j)+p(k-kp+1,i)*p(k,j))*h0(k,ks-kp+1)
  120   continue
      end if
*
      IF (REL) HL=HL-2.d0*RLSHFT(i,j)
      RETURN
  999 continue
      WRITE(IWRITE,4) EL(i),L(i),EL(j),L(j)
4     FORMAT(10X,'UNALLOWED L VALUES OCCURRED IN HL SUBROUTINE'/
     :   2(10X,A3,' HAS L = ',I3))
      stop
      END
      integer function idmax(n,dx,incx)
c
c     finds the index of element having largest value.
c
      double precision dx(1),dmax
      integer i,incx,ix,n
c
      idmax = 0
      if( n .lt. 1 ) return
      idmax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dx(1)
      ix = ix + incx
      do 10 i = 2,n
         if(dx(ix).le.dmax) go to 5
         idmax = i
         dmax = dx(ix)
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dx(1)
      do 30 i = 2,n
         if(dx(i).le.dmax) go to 30
         idmax = i
         dmax = dx(i)
   30 continue
      return
      end
*     ------------------------------------------------------------------
*    		I M P R O V E
*     ------------------------------------------------------------------
*
*     If initial estimates were screened hyrogenics, solve equation
*     without exchange using matrix methods, for lowest eigenfunction
*     satisfying orthogonality.
*
*
      SUBROUTINE improve(ind)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER (nt=100,kx=8,nwd=10)
*
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rkm(nt,nt,kx,kx), lh, krk
*
      LOGICAL fail,omit,rel,all,trace
      COMMON /TEST/fail,omit,rel,all,trace
*
      DIMENSION h1(nt,kx),yrm(nt,kx),a(nt,nt),b(nt,nt),
     :          eivr(nt),eivi(nt),den(nt),eigvec(nt,nt), w1(nt),
     :          eigval(nt), w2(nt),ind(nwd)
      Equivalence (w1(1), eivr(1)), (w2(1),eivi(1)), (eigval(1),den(1))
*
      do 300 i = 1,nwf
	if (ind(i) .ne. -2 ) then
	  node = n(i)-l(i) -1
   	  call potl(i,rel,yrm)
	  call   hlm(z,l(i),nt,kx,ks,ns,db2,rm1,rm2,h1)
	  do 10 m = 1,ks
	    do 10 ii = 1,ns
	      h1(ii,m) = h1(ii,m) - yrm(ii,m)
10        continue
*
*         .. form full coulomb matrix
*
	  call coulom(nt,kx,ks,ns,sb,h1,a,b)
*         .. add orthogonality constraints to inner orbitals
	  icount = 0
	  do 20 j = 1,i-1
	    if( e(i,j) .ne. 0) then
	      icount = icount + 1
 	      call bxv(nt,ks,ns,sb,p(1,j),w1)
 	      call bxv(nt,ks,ns,h1,p(1,j),w2)
	      do 30 jj = 1,ns
	        do 30 ii = 1,ns
	          a(ii,jj) = a(ii,jj) - w1(ii)*w2(jj)
30            continue
	    end if
20        continue
*
*	  if (icount .eq. 0) then
*	    call rsg(nt,ns-2,a(2,2),b(2,2),eigval,1,eigvec,w1,w2)
*	    keig = ns-node-2
*	  else 
	    call rgg(nt,ns-2,a(2,2),b(2,2),eivr,eivi,den,1,eigvec,ier)
	    do 40 ii = 1,ns-2
	      eivr(ii) = eivr(ii)/den(ii)
40          continue
	  print * , 'Eigenvalues:'
	  print '(6f12.8)', (eivr(ii),ii=1,ns-2)
	  print '(6f12.8)', (eivi(ii),ii=1,ns-2)
	  print '(6f12.8)', (den(ii),ii=1,ns-2)
	    node = node -icount
	    do 50 ii = 1,node
	      keig = idmax(ns-2,eivr,1)
	      eivr(keig) = 0.d0
50          continue
	    keig = idmax(ns-2,eivr,1)
	    print *, 'keig, eigval', keig, eivr(keig)
*	  end if
	  do 200 ii = 1,ns-2
	    p(ii+1,i) = eigvec(ii,keig)
200       continue
	  pn=1.d0/sqrt(quadr(i,i,0))
	  if (p(l(i)+3,i) .lt. 0.d0) pn = -pn
	  do 210 ii = 1,ns
	    p(ii,i) = pn*p(ii,i)
210       continue
	  Print *, 'Eigenvector: ',keig
	  print '(6f12.8)', (p(ii,i),ii=1,ns)
	end if
300   continue
      end

*
*     ------------------------------------------------------------------
*               i n i t
*     ------------------------------------------------------------------
*
*       initializes basic constants of the program including those
*   which define the average energy of a configuration.
*
*
      SUBROUTINE init
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER (nt=100,kx=8,nwd=10)
      COMMON /eav/cca(10),ccb(35)
      COMMON /fact/gam(100)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
*  ***** set the factorial needed by rme
*
      call factrl(32)
*
*  *****  set fine structure constant
*
      fine = 0.25d0/(137.036)**2
*
*  *****  set the starting point, step size, and related parameters
*
      rho = -4.d0
      h   = 1./16.d0
      h1 = h/1.5d0
      h3 = h/3.d0
      ch = h*h/12.d0
      eh = dexp(-h)
*
*  *****  average interactions for equivalent electrons
*
*  *****  p - p
*
      cca(1) = 2.d0/25.d0
*
*  *****  d - d
*
      cca(2) = 2.d0/63.d0
      cca(3) = 2.d0/63.d0
*
*  *****  f - f
*
      cca(4) =   4.d0/ 195.d0
      cca(5) =   2.d0/ 143.d0
      cca(6) = 100.d0/5577.d0
*
*  *****  g - g
*
      cca(7) =   20.d0/  1309.d0
      cca(8) =  162.d0/ 17017.d0
      cca(9) =   20.d0/  2431.d0
      cca(10) = 4410.d0/371943.d0
*
*
*  ***** average interactions for non-equivalent electrons
*
*  *****  s - ( s, p, d, f, g )
*
      ccb(1) = 1.d0/ 2.d0
      ccb(2) = 1.d0/ 6.d0
      ccb(3) = 1.d0/10.d0
      ccb(4) = 1.d0/14.d0
      ccb(5) = 1.d0/18.d0
*
*  *****  p - ( p, d, f, g )
*
      ccb(6) = 1.d0/  6.d0
      ccb(7) = 1.d0/ 15.d0
      ccb(8) = 1.d0/ 15.d0
      ccb(9) = 3.d0/ 70.d0
      ccb(10) = 3.d0/ 70.d0
      ccb(11) = 2.d0/ 63.d0
      ccb(12) = 2.d0/ 63.d0
      ccb(13) = 5.d0/198.d0
*
*  *****  d - ( d, f, g )
*
      ccb(14) =  1.d0/ 10.d0
      ccb(15) =  1.d0/ 35.d0
      ccb(16) =  1.d0/ 35.d0
      ccb(17) =  3.d0/ 70.d0
      ccb(18) =  2.d0/105.d0
      ccb(19) =  5.d0/231.d0
      ccb(20) =  1.d0/ 35.d0
      ccb(21) = 10.d0/693.d0
      ccb(22) =  5.d0/286.d0
*
*  *****  f - ( f, g )
*
      ccb(23) =  1.d0/  14.d0
      ccb(24) =  2.d0/ 105.d0
      ccb(25) =  1.d0/  77.d0
      ccb(26) = 50.d0/3003.d0
      ccb(27) =  2.d0/  63.d0
      ccb(28) =  1.d0/  77.d0
      ccb(29) = 10.d0/1001.d0
      ccb(20) = 35.d0/2574.d0
*
*  *****  g - ( g )
*
      ccb(31) =   1.d0/   18.d0
      ccb(32) =  10.d0/  693.d0
      ccb(33) =   9.d0/ 1001.d0
      ccb(34) =  10.d0/ 1287.d0
      ccb(35) = 245.d0/21879.d0
      return
      end
.main

*     ------------------------------------------------------------------
*                       l o o k - t m
*     ------------------------------------------------------------------
*
*     add the deviations to the average energy for a partially filled
*       p- or d- shell
*
      SUBROUTINE looktm( l, sl, sen, q, ip, nsl)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
 
      CHARACTER sl*2, sen*1
      INTEGER iptr(5)
      CHARACTER*3 terms(51)
      data    iptr/6,11,19,35,51/
      data        terms/'3p2','1d2','1s0','4s3','2d3','2p1',
*             .. d2 and d3 terms
     :      '3f2','3p2','1g2','1d2','1s0','4f3','4p3','2h3','2g3',
     :      '2f3','2d1','2d3','2p3',
*            ... d4 terms ...
     :       '5d4','3h4','3g4','3f2','3f4','3d4','3p2','3p4',
     :       '1i4','1g2','1g4','1f4','1d2','1d4','1s0','1s4',
*            ... d5 terms ...
     :       '6s5','4g5','4f3','4d5','4p3','2i5','2h3','2g3',
     :       '2g5','2f3','2f5','2d1','2d3','2d5','2p3','2s5'/
 
*
*  --- search for a partially unfilled p- or d-shell
*
      n = q
      if (n .gt. 2*l+1) n = 4*l+2 - n
      ip = 0
      nsl = 0
      if (n .gt. 1  .and. l .le. 2) then
         if (l .eq. 1) then
            ibegin = 1
            iend = 6
         else
            ibegin = iptr(n-1) + 1
            iend = iptr(n)
         end if
1        i = ibegin
10       if (sl .eq. terms(i)(1:2)) then
            if (sen .eq. ' ' .or. sen .eq. terms(i)(3:3)) then
	       nsl = nsl + 1
	       ip = i
            end if
         end if
         i = i+1
         if (i .le. iend) go to 10
      else if ( n .eq. 1 .and. sl(1:1) .eq. '2') then
	 nsl = 1
      end if
      return
      end
*
*     -----------------------------------------------------------------
*                l o o k - u p
*     -----------------------------------------------------------------
*
      SUBROUTINE lookup(tab,p1,p2,ind,no,key)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      INTEGER tab(*),p1,p2,ind,no,key
      do 40 i = p1,p2
         if (tab(i).eq.key) then
            no = no + 1
            ind = i
         end if
40    continue
      end
*
*     -----------------------------------------------------------------
*                 l v a l
*     -----------------------------------------------------------------
*
*
      INTEGER FUNCTION lval(symbol)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      CHARACTER*1 symbol
      CHARACTER*22 set
      data         set/'spdfghiklmnspdfghiklmn'/
      locate = index(set,symbol)
      if ( locate .le. 11) then
            lval = locate - 1
         else
            lval = locate - 12
      endif
      return
      end
*     ------------------------------------------------------------------
*               m e n u
*     ------------------------------------------------------------------
*
*
*     this routine evaluates a variety of atomic parameters as
*     requested by the user.
*
*
      SUBROUTINE menu
     
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER (nt=100,kx=8,nwd=10)
      INTEGER ouf
      COMMON /inout/ iuf,ouf,iscw
      CHARACTER config*50,el*3,atom*6,term*6
      COMMON /LABEL/config,el(nwd),atom,term
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      LOGICAL fail,omit,rel,all,trace
      COMMON /TEST/fail,omit,rel,all,trace
     
      character*3 el1,el2,el3,el4,func*1,funct*1
     
4     write(iscw,5)
5     format(//5x,
     :  ' these various functions are available:',//
     :  10x,'1 - expectation values of r**k'/
     :  10x,'2 - slater or magnetic integrals'/
     :  10x,'3 - electron density at the nucleus'/
     :  10x,'4 - spin-orbit parameter'/
     :  10x,'5 - transition integrals'/
     :  10x,'6 - exit to main program'/)
      write(iscw,'(5x,a)')'input number 
     :  corresponding to your selection:'
      read(5,'(i1)') ifunc
      go to (10,20,30,40,50,60) ifunc
     
*  ****  compute expectation values
     
10    write(iscw,'(/5x,a,/a,t22,a)')
     :  'input label for electron followed by k: example',
     :  '   2p  3','format(1x,a3,i3) '
      read(5,'(1x,a3,i3)') el1,k
      call eptr(el,el1,i,*10)
      rkev = quadr(i,i,k)
      write(3,12) el1,k,el1,rkev
      write(iscw,12) el1,k,el1,rkev
12    format(/15x,' value of <',a3,'|r**',i2,'|',a3,'> = ',1pd14.7,
     :         ' a.u.'/)
      go to 4
     
*  ****  determine slater integrals  fk, gk, rk, nk, mk, vk
     
     
20    write(iscw,'(/5x,a/a,t22,a)')
     : 'input parameters for  fk,gk,rk,nk,mk or vk integral: example',
     : ' f 0( 1s, 2s)','format:  (a1,i2,1x,4(a3,1x)) '
      read(5,'(a1,i2,1x,4(a3,1x))') funct,k,el1,el2,el3,el4
      if (funct .ge. 'a' .and. funct .le. 'z')
     :   func = char(ichar(funct) + ichar('a') - ichar('a'))
      call eptr(el,el1,i1,*20)
      call eptr(el,el2,i2,*20)
      if ( el3 .ne. ' ') call eptr(el,el3,i3,*20)
      if ( el4 .ne. ' ') call eptr(el,el4,i4,*20)
      if (func .eq. 'f' ) then
        si = fky(i1,i2,k,rel)
      else if (func .eq. 'g') then
        si = gky(i1,i2,k,rel)
      else if (func .eq. 'r') then
        si = rky(i1,i2,i3,i4,k,rel)
      else if (func .eq. 'n') then
*       si = sn(i1,i2,i2,i1,k)
      else if (func .eq. 'm') then
*       si = sn(i1,i2,i1,i2,k)
      else if (func .eq. 'v') then
*       si = vk(i1,i2,i2,i1,k) - vk(i2,i1,i1,i2,k)
      else
        write(iscw,41)
41      format(15x,'integral unknown: re-enter')
        go to 20
      end if
      if (func .ne. 'r') then
         write(3,25) func,k,el1,el2,si,219474.d0*si
         write(iscw,25) func,k,el1,el2,si,219474.d0*si
25       format(/15x,
     :  'integral  ',a1,i2,'(',a3,',',a3,') = ',1pd14.7,' a.u.'/
     :   40x,0pf14.3,' cm-1'/)
      else
         write(3,26) func,k,el1,el2,el3,el4,si,219474.d0*si
         write(iscw,26) func,k,el1,el2,el3,el4,si,219474.d0*si
26       format(/15x,
     :  'integral  ',a1,i2,'(',2a3,',',2a3,') = ',1pd14.7,' a.u.'/
     :   46x,0pf14.3,' cm-1'/)
      end if
      go to 4
     
*  ****  compute electron density at the nucleus
     
30    write(iscw,'(/5x,a/a,t22,a)')
     :  'input identifying label for electron: example',
     :  '   1s','format(1x,a3) '
      read(5,'(1x,a3)') el1
      call eptr (el,el1,i,*30)
      ll = l(i)
      if (ll .eq. 0) then
         d = az(z,h,ks,l(i)+1,p(l(i)+2,i))**2
      else
         d = 0
      end if
      write(3,32) el1,d
      write(iscw,32) el1,d
32    format(/15x,'density at the nucleus for ',a3,' = ',1pd14.7,
     :     ' a.u.'/)
      go to 4
     
*  ****  compute spin-orbit parameter
     
40    write(iscw,'(/,5x,a/a,t22,a)')
     :  'input identifying label for electron: example',
     :  '   2p','format(1x,a3) '
      read(5,'(1x,a3)') el1
      call eptr (el,el1,i,*40)
      zeta = 0.d0
      if (l(i) .ne. 0)  zeta = bwzeta(i)
      zetacm = 219474*zeta
      write(3,43) el1,zeta,zetacm
      write(iscw,43) el1,zeta,zetacm
43    format(/15x,'spin-orbit parameter for ',a3,' = ',1pd14.7,
     :        ' a.u.'/46x,0pf14.3,' cm-1'/)
      go to 4
     
*  ****  compute transition integrals
     
50    write(iscw,'(/5x,a/a,t22,a)')
     :  'input identifying labels and power of r: example',
     :  ' t 1( 2s, 2p)','format:  (a1,i2,2(1x,a3)) '
      read(5,'(a1,i2,1x,a3,1x,a3)') func,k,el1,el2
      call eptr (el,el1,i1,*50)
      call eptr (el,el2,i2,*50)
      ti = quadr(i1,i2,k)
      write(3,52) func,k,el1,el2,ti
      write(iscw,52) func,k,el1,el2,ti
52    format(/15x,
     :  'integral  ',a1,i2,'(',a3,',',a3,') = ',1pd14.7,' a.u.'/)
      go to 4
60    return
      end
*     ------------------------------------------------------------------
*               m e t h o d
*     ------------------------------------------------------------------
*
*       uses m1, m2, or m3 to solve the radial equation. if the input
*   data indicated meth(i) = 3, then this  solution  is  returned  to
*   de.  otherwise,  the routine searches for an acceptable  solution
*   which  is  both  positive  near  the  origin and has the required
*   number  of nodes.
*
      SUBROUTINE method(i,first,ed,pde,m)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(nt=100,kx=8,nwd=10)
*
      COMMON /coeff/coef(100),ijptr(5,5)
      character config*50,el*3,atom*6,term*6
      COMMON /label/config,el(nwd),atom,term
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/ sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :             iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      LOGICAL fail,omit,rel,all,trace
      COMMON /test/fail,omit,rel,all,trace
      DOUBLE PRECISION pde(nt),h1(nt,kx),bpv(nt),w(nt),yrm(nt,kx),xv(nt)
      LOGICAL first
      SAVE em,eu,fu,fm
*
*  *****  'first' must be 'true' the first time solve is called for
*  *****  potential and exchange to be computed
*  *****  'eu' is the upper bound of the energy parameter
*  *****  'em' is the minimum value of the energy parameter
*
      fail = .false.
      mit = 0
      if (first) then
        first = .false.
        em = 0.d0
        eu = ((z - dmin1(d5*s(i),d2*s(i)))/n(i))**2
        fu = eu
 	call potl(iel,rel,yrm)
	call xch(iel,rel,xv)
	call hlm(z,l(iel),nt,kx,ks,ns,db2,rm1,rm2,h1)
	do 10 m = 1,ks
	  do 10 ii = 1,ns
	    h1(ii,m) = h1(ii,m) - yrm(ii,m)
10      continue
	call bxv(nt,ks,ns,h1,p(1,iel),bpv)
	do 12 ii = 1,ns
	  bpv(ii) = bpv(ii) - xv(ii)
12      continue
	ed = ddot(ns,p(1,iel),1,bpv,1)
	do 14 j = 1,nwf
	  eij = e(iel,j)
	  if( eij .ne. 0 .and. iej .ne. 1.d-10 .and. iel.ne.j) then
	    call bxvpw(nt,ks,nc,eij,sb,p(1,j),xv,w)
	  end if
14	continue
      end if
17    call solve(i,first,ed,h1,xv,pde,m)
*
*  *****  if mk equals 3, omit the node checking
*
      if (mk .le. 2) then
*
*  *****  count the number of nodes
*
        mn = m
        nc = nodec(mn)
        if (trace) write (6,99) el(i),nc,mn,pde(mn),ed,eu,em
99      format(2x,a3,' nc =',i3,' mn =',i3,' pde(mn) =',
     :         d10.2,' ed =',d10.2,' eu =',d10.2,' em =',d10.2)
*
*  *****  if node count is off by no more than 1 and deltae is still
*  *****  quite large, apply the deltae correction
*
        if (iabs(nc-node) .eq. 0) then
*
*         .. we have the correct number of nodes
          if (pde(mn) .lt. 0.d0 ) then
            do 52 j = 1,no
              pde(j) = - pde(j)
52          continue
	  end if
          return
	else if (nc - node .lt. 0) then
*
*  *****  the solution has too few nodes
*
          if (pde(mn) .gt. 0.d0) then
            del = 1.d0 - ed/eu
            eu = ed
            if ( del .lt. .05d0) fu = fu*((l(i)+1+nc)/fn)**2.5
            if (del  .ge. .05d0) fu = ed*((l(i)+1+nc)/fn)**2.5
            if (fu .lt. em) fu = d5*(eu + em)
            if (dabs(fu - ed) .lt. 0.001d0) then
              write (6,28) mk,el(i),nc,ed,em,eu
	      fail = .true.
	      return
	    else
	      mit = mit+1
              ed = fu
            end if
	  else
*
*           ..try a new value of ed which must lie within the upper and
*             lower bound
*
            edp = ed
            ed = ed*((l(i)+1+nc)/fn)**2.5
            if (ed .ge. eu ) ed = d5*(eu + edp)
            if (ed .le. em ) ed = d5*(em + edp)
33          mit = mit + 1
	  end if
	else 
*
*         ..the solution has too many nodes
*
          if (pde(mn) .gt. 0.d0) then
            del = 1.d0 - em/ed
            em = ed
            if (del .lt. 0.05d0) fm = fm*((l(i)+1+nc)/fn)**2.5
            if (del .ge. 0.05d0) fm = ed*((l(i)+1+nc)/fn)**2.5
            if (fm .gt. eu) fm = d5*(eu + em)
            if (dabs(fm - ed) .lt. 0.001d0) then
               write (6,28) mk,el(i),nc,ed,em,eu
	       fail = .true.
	       return
	    end if
            ed = fm
          end if
	end if
      end if
      if ( eu .le. em ) write(6,30) em,eu,ed
30    format(6x,48hwarning: difficulty with node counting procedure/
     :       6x,42hlower bound on ed greater than upper bound/
     :       6x,5hel = ,f10.6,7h  eu = ,f10.6,7h  ed = ,f10.6)
      if ( mit .gt. 3*n(i) .or. eu-em .lt. fn**(-3)) then
         write (6,28) mk,el(i),nc,ed,em,eu
	 fail = .true.
	 return
      else
        go to 17
      end if
*
28    format(10x,6hmethod,i2,38h unable to solve equation for electron,
     :   a3/10x,5hnc = ,i3,3x,5hed = ,f10.6,3x,5hel = ,
     :   f10.6,3x,5heu = ,f10.6)
      end
*
*     ------------------------------------------------------------------
*               n o d e c
*     ------------------------------------------------------------------
*
*      counts the number of nodes of the function pde(j) in the range
*   j = 40,...,m-10.   the node counting procedure counts the local max
*   and min values.   only nodes between sufficiently large max and
*   min values are counted.
*
*
      INTEGER FUNCTION nodec(pde,ns,m)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION pde(*)
*
*   ***** find max |pde(j)|
*
      m = idamax(ns,pde,1)
      pdem = pde(m)
*
*   *****  count the number of local max or min's
*
      ncc = 0
      sign = 0.d0
      diff1 = pde(4) - pde(3)
      do 2 j = 4,ns-1
        diff2 = pde(j+1) - pde(j)
        if (diff2*diff1 .lt. 0.d0 .and. diff1 .ne. 0.d0) then
*
*   *****  a max or min has been found.   test if it is
*          sufficiently large
*
          if ( dabs(pde(j))/pdem .gt. 0.05d0 ) then
*
            if (sign .eq. 0.d0 ) then
*              .. this is the first significant maximum
               m = j
	    else
*           .. this is not the first, test whether a sign change has
*              occurred since the last significant max or min
              if (pde(j)*sign .lt. 0.d0 ) ncc = ncc + 1
	    end if
*           .. reset for the next anti-node
            sign = pde(j)
          end if
	end if
        diff1 = diff2
2     continue
      nodec = ncc
      end
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
*
*     ------------------------------------------------------------------
*               o u t p u t
*     ------------------------------------------------------------------
*
*       the radial functions and orthogonality integrals are printed,
*   if print is .true.   the  functions  will  also  be  punched  (or
*   stored) on unit ouf, if ouf .ne. 0.
*
      SUBROUTINE output(print)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(nt=100,kx=8,nwd=10)
      LOGICAL print
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
      CHARACTER config*50,el*3,atom*6,term*6
      COMMON /LABEL/config,el(nwd),atom,term
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      dimension out(8)
      if ( print ) then
*
*  *****  print radial functions, 7 per page
*
        ml = ib
2       mu = min0(ml+7,nwf)
        i = mu - ml + 1
        mx = 0
        do 1 j = ml,mu
1         mx = max0(mx,max(j))
        write (3,5) atom,term,(el(j),j=ml,mu)
5       format(1h1,9x,19hwave functions for  ,2a6//10x,1hr,8(10x,a3))
        k= 0
        kk = 0
        do 6 j = 1,mx
          do 9 jj = ml,mu
            ij = jj - ml + 1
            out(ij) = p(j,jj)
 9        continue
          k = k+1
          if (k .lt. 10) then
            k = 1
            kk = kk+1
            if (kk .ge. 5) then
              kk = 0
              write (3,23)
23            format(1h1//)
            else
              write (3,10) t(j+ks),(out(jj),jj=1,i)
10            format(/(f13.5,f15.6,7f13.6))
	    end if
          end if
 6      continue
        do 15 j = ml,mu
          ij = j - ml + 1
          out(ij) = dpm(j)
15      continue
        write (3,16) (out(j),j=1,i)
16      format(4x,10hmax. diff. ,f15.7,7f13.7)
        ml = ml+8
        if (ml .le. nwf) go to 2
      end if
30    if ( ouf .ne. 0) then
*
*  *****  output functions on unit ouf for future input
*
        do 3 i = 1,nwf
          mx = max(i)
          write (ouf) atom,term,el(i),mx,z,e(i,i),
     :                (p(j,i),j=1,mx)
3       continue
      end if
*
      end
*     ------------------------------------------------------------------
*               p o t l
*     ------------------------------------------------------------------
*
*       computes and stores the potential function
*                              2(k-1)
*              yr = sum  a    y      (j,j;r)
*                   j,k   ijk
*
      SUBROUTINE potl(i,rel,yrm)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(nt=100,kx=8,ktx=3*kx-2,nwd=10)
      DOUBLE PRECISION yrm(nt,kx)
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rkm(nt,nt,kx,kx),lh,krk

      COMMON /SLATER/fc1(nt,kx),fc2(nt,kx),yk(nt),fyk(nt,kx),
     :        dyk(ktx,nt),ipvtd(nt),fc(nt,kx),ic1,ic2,iy1,iy2,ky
*
      LOGICAL kfirst, REL
      DOUBLE PRECISION yr(nt,kx), yv(nt)
*---------------------------------------------------------------------
        print *, 'Entering Potl for orbital ',i
        nint = ns-ks+1
*       .. clear yv
        do 10 ii = 1,ns
	  yv(ii) = 0.d0
   10   continue
*       .. clear yr
	do 15 m = 1,ks
	  do 15 ii = 1,ns
	    yr(ii,m) = 0.d0
   15   continue
*
*       .. compute contribution from each k
        do 20  k = 0,2*lmax,2
	  ykn   = 0.0d0
	  kfirst = .true.
	  do 30 j = 1,nwf
	    if (k .le. 2*min(l(i),l(j))) then
              cc = 2.d0*a(i,j,k)
	      print *, 'i,j,k,a(i,j,k)',i,j,k,cc
	      if (abs(cc) .gt. 1.d-10) then
		if (k .eq. 0) ykn = ykn + cc
		if ( kfirst ) then
		   kfirst = .false.
*                  .. clear array for fyk
		   do 40 m = 1,ks
		     do 40 jj = 1,nint
		       fyk(jj,m) = 0.d0
 40		   continue
*                  .. clear the array yk
	           do 45 ii=1,ns
	             yk(ii)=0.d0
 45 	           continue
*
*                  set up and factor the differential operator
*
		   print *, ' Potl: calling facdyk, k=',k
                   call facdyk(nt,kx,ktx,ks,k,nint,t,bsp,bspd,
     :                         db2,rm2,ipvtd,dyk)
		end if
*   
	        if (ic1 .ne. j) then
	          call yval(0,0,nt,kx,ks,nint,gr,bsp,p(1,j),fc1)
	          ic1 = j
	        end if
		do 50 m = 1,ks
		  do 50 ii = 1,nint
		    fyk(ii,m) = fyk(ii,m) + cc*fc1(ii,m)*fc1(ii,m)
 50             continue
	      end if
	    end if
 30       continue
*
*         .. if there is a contribution, add
*
          if (.not. kfirst) then
*
*         add a relativistic correction if rel=TRUE
*
	    if (rel) then
	      c = (2*k+1)*fine
              do 60 m = 1,ks
	        do 60 ii = 1,nint
	          yr(ii,m)= yr(ii,m)+ c*fyk(ii,m)*grm(ii,m)
  60          continue
	    end if
*
*   set up the array yk (spline approx. of SUM ()(1/r)fc1(r)fc2(r)
*
	    do 70 m=1,ks
	      do 70 ith=1,ks
	        do 70 ii = 1,nint
     	        irow=ii+ith-1
	          yk(irow)=yk(irow)+grw(ii,m)*grm(ii,m)*bsp(ii,m,ith)*
     :  	        fyk(ii,m)
   70	    continue
*
    	    const=-(2*k+1.d0)
	    do 80 ii=1,ns
	      yk(ii)=const*yk(ii)
   80	    continue
	    yk(1)=0.d0
	    yk(ns) = ykn
*           print *, ' Potl: yk'
*           print '(6f12.8)', (yk(ii),ii=1,ns)
*
*   solve the matrix equation
*
	    call dgbsl(dyk,ktx,ns,ks-1,ks-1,ipvtd,yk,0)
*           print *, 'Potl: after dgbsl', i
*           print '(6f12.8)', (yk(ii),ii=1,ns)
*
*            .. add to yv
	    do 90 ii = 1,ns
	      yv(ii) = yv(ii) + yk(ii)
 90         continue
	  end if
 20     continue 
*      print *, 'Potl: yv', i
*      print '(6f12.8)', (yv(ii),ii=1,ns)
*
*  evaluates the function Y(r) at all the gaussian points
*
        call yval(1,0,nt,kx,ks,nint,grm,bsp,yv,yr)
	print *, ' Y(r)'
	do 100 m = 1,ks
	  do 100 ii = 1,nint
 	    yr(ii,m) = yr(ii,m)*grw(ii,m)*grm(ii,m)
 100    continue
	call minty(0,nt,kx,ks,nint,bsp,yr,yrm)
        print *, ' Test yrm'
	call bxv(nt,ks,ns,yrm,p(1,i),yv)
        ans1 = ddot(ns,yv,1,p(1,i),1)
	print *, '<p1,potl,p1>', ans1, quadr(1,1,-1)
*
	end
*
*     ------------------------------------------------------------------
*              Q U A D
*     ------------------------------------------------------------------
*
*       Returns the value of <a, S b> where a and b are spline
*     expansion coefficients
*
      DOUBLE PRECISION  FUNCTION quad(ns,ks,a,b)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER (nt=100,kx=8,nwd=10)
      DOUBLE  PRECISION a(nt), b(nt)
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rk(nt,nt,kx,kx),lh,krk
*
      quad = 0.d0
      do 100 i = 2,ns
	quad = quad + a(i)*b(i)*sb(i,ks)
  100 continue
      do 120 m = 2,ks
*       ..i-m+1  < m
	do 120 i = m+1,ns
	  quad = quad +(a(i)*b(i-m+1)+a(i-m+1)*b(i))*sb(i,ks-m+1)
  120 continue
      print *, 'quad', quad
*
      END
*=======================================================================
        subroutine refine(nt,n,a,b,c,eigval,eigvec,ipvt,x,ierr)
*=======================================================================
*   Given an estimate of the eigenvalue, 
*           (a - lambda b)v = c v = 0
*   computes and eigenvector and updates the eigenvalue
*
*   Written by C. F. Fischer
*------------------------------------------------------------------------
*
	implicit double precision (a-h,o-z)
	dimension a(nt,nt), b(nt,nt), c(nt,nt),
     :            eigvec(nt),ipvt(nt),x(nt)
*
*   on entry
*   --------
*     nt     dimensions of a, b
*     n       size of the matrix
*     a,b     matrices
*     eigval  estimate of lambda
*     ipvtd   integer work array
*     x       double precision work array
*
*   on exit
*   -------
*     eigval  improved estimate of the eigenvalue
*     eigvec  the eigenvector corresponding to eigval
*
*-----------------------------------------------------------------------
*
* 	..Form the c matrix
	s = eigval
	do 10 j = 1,n
	   do 10 i = 1,n
	     c(i,j) = a(i,j) -s*b(i,j)
   10   continue
*
*       .. form an initial guess of the eigenvector
	do 20 i = 1,n
	  eigvec(i) = 0.d0
   20   continue
	eigvec(2) = .5d0
	eigvec(3) = .5d0
	eigvec(4) = .5d0
	eigvec(5) = 1.d0

   30   call ginvit(c,nt,n,eigvec,ipvt,x,ierr)
*
	eiv = 0.d0
	den = 0.d0
	do 50 i = 1,n
	  eiv = eiv + a(i,i)*eigvec(i)**2
	  den = den + b(i,i)*eigvec(i)**2
   50   continue
*
*       .. the off-diagonal part
	eiv2 = 0.d0
	den2 = 0.d0
	do 60 j = 1,n-1
	  t = 0.d0
	  d = 0.d0
	  do 70 i = j+1,n
             t= t + a(i,j)*eigvec(i)
	     d= d + b(i,j)*eigvec(i)
   70     continue
	  eiv2 = eiv2 + t*eigvec(j)
	  den2 = den2 + d*eigvec(j)
   60   continue
 
	eivnew = (eiv + 2*eiv2)/(den +2*den2)
	eigval = eivnew
*
*       .. normalize the eigenvector
	dn = 1.d0/sqrt(den+2*den2)
	do 80 i = 1,n
	  eigvec(i) = eigvec(i)*dn
   80   continue

	return
	end
*
*     ------------------------------------------------------------------
*               r e f o r m
*     ------------------------------------------------------------------
*
*     convert the free-format str1 to the fixed 5(1x,a3,1x,i4,1x) format
*     for str2
*
      SUBROUTINE reform(str1,str2)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
      CHARACTER*50 str1,str2,blank
      data                   blank/'   '/
*
    1 i = 0
      str2 = blank
      is = 0
    2 js = index(str1(is+1:),'(')
      if (js .ne. 0) then
         if (js .gt. 5) go to 10
         i = i+5
         str2(i-js+1:i) = str1(is+1:is+js)
         is = is + js
         js = index(str1(is+1:),')')
         if (js .eq. 0 .or. js .gt. 5) go to 10
         i = i+5
         str2(i-js+1:i) = str1(is+1:is+js)
         is = is + js
         go to 2
      end if
      return
   10 write(iscw,*)' error in ',str1,': re-enter'
      read(5,'(a)') str1
      go to 1
      end
*
*   --------------------------------------------------------------------
*               r e o r d
*   --------------------------------------------------------------------
*
*       reorder the list of first appearance so that the FUNCTIONs to be
*   iterated appear last in the list.
*
        SUBROUTINE reord(of, elc, nwf, ierr)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
        CHARACTER*3 of(30), elc
*
        ierr = 1
        call eptr(of, elc, i, *99)
        do 10 j = i, nwf-1
           of(j) = of(j+1)
10      continue
        of(nwf) = elc
        ierr = 0
99      return
        end
*     ==================================================================
*              R E L - S H I F T
*     ==================================================================
*
*
*       Computes the mass velocity  and one-body Darwin term
*   corrections for the relativistic shift in the energy of the electron
*   including non-diagonal corrections
*
*
      DOUBLE PRECISION  FUNCTION RLSHFT(I1,I2)
      IMPLICIT DOUBLE PRECISION (a-h,o,p,r-z)
      PARAMETER (nt=100,kx=8,nwd=10)
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rk(nt,nt,kx,kx),lh,krk
*
      DIMENSION vc(nt,kx),a(nt,kx)
*
* *****  form the mass-velocity correction matrix in the spline basis
*
      LL = L(i1)
      print *, i1,i2,ll
      call mvc(LL,nt,kx,ks,ns,grm,grw,bsp,bspd,vc)
*
* .... store the product of coefficients
*
      last = min0(max(i1),max(i2))
      print *, last
      do 10 j =1,last
        a(j,1) = p(j,i1)*p(j,i2)
  10  continue
      do 20 jp = 2,ks
	do 20 j = 1,last-jp+1
	  a(j,jp) = p(j,i1)*p(j+jp-1,i2) + p(j,i2)*p(j+jp-1,i1)
  20  continue
*
* ...  form the sums
*
	r = 0.d0
        do 100 j = 2,last
	  r = r + a(j,1)*vc(j,ks)
  100   continue
	do 110 jp = 2,ks
*         .. j' < j
	  do 110 j = jp+1,last
	    r = r + a(j-jp+1,jp)*vc(j,ks-jp+1)
  110     continue

      print *, 'RLSHFT: relsh,z,h,ks,last',r,z,h,ks,last,fine
      IF ( L(I1) .EQ. 0 ) then
	RELSH = -r + Z*AZ(z,h,ks,LL+1,p(2,I1))*AZ(z,h,ks,LL+1,p(2,I2))
      else
	RELSH = -r
      end if
      print *, 'RLSHFT:', relsh
      RLSHFT = 0.5d0*RELSH*FINE
      RETURN
      END
*     ------------------------------------------------------------------
*                r m e
*     ------------------------------------------------------------------
*
*
      double precision FUNCTION rme(l,lp,k)
*
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
*
      COMMON/fact/gam(100)
*
*--- evaluates the reduced matrix element (l//c(k)//lp)  -  see fano
*    and racah, irreducible tensorial sets, chap. 14, p. 81
*
      if (min0(l,lp) .eq. 0) then
         rme = 1.d0
       else if ( k .eq. 0) then
         rme = 2*l+1
         rme = dsqrt(rme)
       else
         i2g=l+lp+k
         ig=i2g/2
         if (i2g - 2*ig .ne. 0) then
             rme = 0.d0
         else
            i1=ig-l
            i2=ig-lp
            i3=ig-k
            qusqrt=(2*l+1)*(2*lp+1)
            rme=dsqrt(qusqrt)*dexp((gam(2*i1+1)+gam(2*i2+1)+gam(2*i3+1)-
     :        gam(i2g+2))/2.d0 +gam(ig+1)-gam(i1+1)-gam(i2+1)-gam(i3+1))
         end if
      end if
      return
      end
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
*
*     ------------------------------------------------------------------
*               s o l v e
*     ------------------------------------------------------------------
*
*       when first is .true., solve computes the potential and exchange
*   function and initializes variables for the i'th radial  equation.
*   the vector p1 is the solution of the radial equation and p2 the
*   variation of the solution with respect to the energy parameter
*   e(i,i).
*
*
      SUBROUTINE solve(iel,first,ed,h1,xv,pde,m)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(nt=100,kx=8,ktx=3*kx-2,nwd=10)
      character config*50,el*3,atom*6,term*6
      COMMON /label/config,el(nwd),atom,term
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/ sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :             iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      LOGICAL fail,omit,rel,all,trace,print,strong
      COMMON /TEST/fail,omit,rel,all,trace
      LOGICAL first
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rkm(nt,nt,kx,kx), lh, krk
*
       DOUBLE PRECISION bandde(ktx,nt),h1(nt,kx),xv(nt),w(nt),ipvtd(nt),
     :                  pde(nt)
*--------------------------------------------------------------------
*
*  clear the bandde array
*
  	do 20 i=1,ks-1
	  do 20 j=1,ns
            bandde(i,j)=0.d0
 20	continue 
*
*  set up bandde
*
*       ... lower portion
	do 30 j=1,ks
   	  do 30 i=ks-j+1,ns
	    bandde(3*ks-1-j,i-ks+j)=h1(i,j) -ed*sb(i,j)
 30	continue            
*      ... upper portion
	do 40 j=2,ks
	  do 40 i=1,ns-j+1
	    bandde(2*ks-j,i+j-1)= h1(i+j-1,ks-j+1)-ed*sb(i+j-1,ks-j+1)
 40	continue
*
* apply the zero boundary condition at the origin and at rmax
*
	do 50 i=1,ks
*           .. diagonals
      	    j=2*ks-i
	    bandde(j,i)=0.d0
	    j=3*ks-1-i
	    bandde(j,ns-ks+i)=0.d0
*           .. columns
	    bandde(2*ks -2 +i,1) = 0.d0
	    bandde(ks-1+i,ns) = 0.d0
 50	continue
 	bandde(2*ks-1,1)=1.d0 
        bandde(2*ks-1,ns) = 1.d0
     
*	print *, ' Band DE'
*	do 99 m = 1,3*ks-2
*	  print '(6f12.8)', (bandde(m,ii),ii=1,ns)
*99      continue
*	print *, ' h1'
*	do 199 m = 1,ks
*	  print '(6f12.8)', (h1(ii,m),ii=1,ns)
*199     continue
*
* factorize bandd
*
	call dgbfa(bandde,ktx,ns,ks-1,ks-1,ipvtd,ier)
	if (ier .ne. 0) then
          print *, 'bandDE is non_invertible'
	  print *, ' Ier = ',ier
	  stop
	end if
*       
*       copy xv to pde
	do 60 i = 1,ns
	  pde(i) = xv(i)
60      continue
	call dgbsl(bandDE,ktx,ns,ks-1,ks-1,ipvtd,pde,0)
* 
	if (meth(i) .eq. 2) then
	  do 70 i = 1,ns
	    w(i) = p(i,iel)
70        continue
	  call dgbsl(bandDE,ktx,ns,ks-1,ks-1,ipvtd,w,0)
	  a00 = quad(ns,pde,pde)
	  a01 = quad(ns,pde,w)
	  a11 = quad(ns,w,w)
	  b = 2.d0*a01/a11
          c = (a00 - 1.d0)/a11
	  if (b .gt. 0.d0) then
	    den = b + sqrt(b*b-4*c)
	  else
	    den = b - sqrt(b*b-4*c)
          end if
	  deltae = - 2.d0*c/den
	  do 80 i = 1,ns
	    pde(i) = pde(i) + deltae*w(i)
80        continue
        end if
	m = ns
	end	
*     ------------------------------------------------------------------
*               s u m m r y
*     ------------------------------------------------------------------
*
*       the results of a calculation are summarized.   these include
*   the following for each electron:
*
*          e(nl)   - diagonal energy parameter
*          i(nl)   - -(1/2)<nl|l|nl>
*          ke      - i(nl) + z <r>
*          rel     - relativistic shift (mass-velocity, darwin term,
*                    spin-spin contact term)
*          sigma   - screening parameter as defined by eq. (6-  ).
*          az(nl)  - starting parameter, p(r)/r**(l+1) as r -> 0.
*          1/r**3  - expected value of <1/r**3>
*          1/r     - expected value of <1/r>
*          r       - expected mean radius
*          r**2    - expected value of <r**2>
*
*   these results are followed by:
*
*          kinetic energy (ek)
*          potential energy (ep) = et - en
*          ratio                 = ep/en
*          non- relativistic energy (et - erel)
*          relativistic shift (erel) for the state
*          total energy (et)
*
      SUBROUTINE summry(et,erel)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(nt=100,kx=8,nwd=10)
      CHARACTER config*50,el*3,atom*6,term*6
      COMMON /coeff/coef(100),ijptr(5,5)
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
      COMMON /LABEL/config,el(nwd),atom,term
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      LOGICAL fail,omit,rel,all,trace
      COMMON /TEST/fail,omit,rel,all,trace
      DIMENSION r1(nwd),rm1(nwd),ek(nwd)
*
      pi = acos(-1.d0)
      write (3,9) atom,term
9     format(/// 24x,'atom ',a6,3x,'term ',a6//
     :   2x,'nl',7x,'e(nl)',
     :   8x,'i(nl)',7x,'ke(nl)',8x,'rel(nl)',3x,'s(nl)',7x,'az(nl)')
      en = 0.d0
      rel = .false.
*
*  *****  compute and print one-electron parameters
*
      do 10 i = 1,nwf
        r1(i) = quadr(i,i,1)
        ek(i) = -0.5d0*hl(el,i,i,rel)
        rm1(i) = quadr(i,i,-1)
        ekinp = ek(i) + z*rm1(i)
        en = en+ sum(i)*ekinp
	print *, i,sum(i),ek(i),rm1(i),ekinp
        rh = 3*n(i)*n(i) - l(i)*(l(i) + 1)
        sc = z - 0.5d0*rh/r1(i)
        s(i) = sc
        rels = rlshft(i,i)
        azz = az(z,h,ks,l(i)+1,p(l(i)+2,i))
        write (3,15)el(i),e(i,i),ek(i),ekinp,rels,s(i),azz
15      format(1x,a3,f14.7,3f13.6,f8.3,f14.6)
10    continue
*
*  *****  compute moments
*
      write(3,8) 'delta(r)'
 8    format(//2x,'nl',6x,a8,5x,'1/r**3',7x,'1/r',9x,'r',8x,'r**2')
      do 11 i = 1,nwf
        rm3 = 0
        if (l(i) .ne. 0) rm3 = quadr(i,i,-3)
        rp2 = quadr(i,i,2)
        rz = 0.
        if ( l(i) .eq. 0) rz = az(z,h,ks,l(i)+1,p(l(i)+2,i))**2/(4.*pi)
        write(3,16) el(i),rz,rm3,rm1(i),r1(i),rp2
16      format(1x,a3,f14.3,f13.4,f11.5,f10.5,f11.5)
11    continue
31    etn = et - erel
      epotl = etn - en
      ratio = epotl/en
      write(iscw,26) etn,en,erel,epotl,et,ratio
      write(3,26) etn,en,erel,epotl,et,ratio
26    format(//5x,'total energy (a.u.)'/5x,'----- ------'/
     : 10x,' non-relativistic   ',f15.8,t50,'kinetic   ',f15.8/
     : 10x,' relativistic shift ',f15.8,t50,'potential ',f15.8/
     : 10x,' relativistic       ',f15.8,t50,'ratio     ',f15.9)
13    return
      end
*     ==================================================================
*    		W A V E F N
*     ==================================================================
*
*       This routine initializes radial functions by the procedure
*   indicated by IND(I).
*
*         Value of IND(I)     Method
*         ---------------     ------
*             -1           Functions read from unit IU2
*              0           Screened hydrogenic functions with ZZ=Z-S(I)
*              1           Functions in memory left unchanged
*                                                  0
*   The set of functions are then orthogonalized, Y (i, i;r) and the
*   diagonal energy parameters computed, when necessary.
*
*
      SUBROUTINE WAVEFN(ind)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER (nt=100,kx=8,nwd=10)
*
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
*
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
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rkm(nt,nt,kx,kx), lh, krk
*
      LOGICAL fail,omit,rel,all,trace,print,strong
      COMMON /TEST/fail,omit,rel,all,trace
*
      DIMENSION ZZ(NWD),PT(NT),ind(nwd),w3(nt),bs(kx,nt),
     :          w21(nt,kx),w22(nt,kx)
      CHARACTER*3 EL1
      CHARACTER*6 AT,TT,ATM(NWD),TRM(NWD)
      CHARACTER*24 TITLE
      LOGICAL LD
*
      nint=ns-ks+1
      call facsb(nt,kx,ks,ns,sb,bs)
*
*  ***** READ THE WAVEFUNCTIONS
*
      IF (IUF .EQ. 0) GO TO 5
2     READ(IUF,END=5) AT,TT,EL1,MM,ZT,ETI,(PT(J),J=1,MM)
      M = min(ns,mm)
      CALL EPTR(EL,EL1,I,*2)
      IF ( I .GT. 0 .AND. IND(I) .EQ. -1) THEN
         ATM(I) = AT
         TRM(I) = TT
         MAX(I) = M
         ZZ(I)  = ZT
         C = 1.d0
         IF ( Z .NE. ZT ) C = Z/ZT
*
*  *****  SCALE RESULTS IF DATA IS FOR AN ATOM WITH A DIFFERENT Z
*
         e(I,I) = C*C*ETI
         DO 11 J = 1,M
            P(J,I) = C*PT(J)
11       CONTINUE
*
*  *****  SET REMAINING VALUES IN THE RANGE = 0.
*
         IF ( M .EQ. ns ) GO TO 12
         M = M +1
         DO 13  J=M,ns
13       P(J,I) = 0.d0
12       IND(I) = -2
      ENDIF
      GO TO 2
*
*  *****  SET PARAMTERS FOR ELECTRONS AND INITIALIZE FUNCTIONS
*
5     continue
      DO 9 I = 1,NWF
      IF (IND(I)) 7,8,9
*
*  ***** WAVE FUNCTIONS NOT FOUND IN THE INPUT DATA, SET IND = 0
*
7     IF ( IND(I) .EQ. -2 ) GO TO 4
      IND(I) = 0
      WRITE(iscw,27) EL(I)
27    FORMAT(8X,'WAVE FUNCTIONS NOT FOUND FOR ',A3)
*
*  *****  DETERMINE ESTIMATES OF THE WAVE FUNCTIONS BY THE SCREENED
*  *****  HYDROGENIC APPROXIMATION
*
8     continue
      CALL BHWF(N(I),L(I),Z-s(i),nt,kx,ks,nint,gr,grw,bsp,bs,w21,w22,
     :          P(1,I))
      M = NS-1
30    IF ( DABS(P(M,I)) .LT. 1.D-15 ) then
        P(M,I) = 0.d0
        M = M-1
        GO TO 30
      END IF
*31    MAX(I) = M+1
 31    MAX(I) = ns
       print *, ' nl expansion', el(i)
       print '(6f12.8)', (p(ii,i),ii=1,ns)
      e(i,i) = ((Z-s(i))/n(i))**2
       print *, 'After orthogonalizationn'
*
*  *****  ORTHOGONALIZE TO INNER FUNCTIONS
*
4     IM = I - 1
      m = max(i)
      DO 6 II =1,IM
	if (e(i,ii) .ne. 0.d0) then
          PN = QUADR(I,II,0)
          print *, 'Overlap between ',i,ii,pn,e(i,ii)
          IF ( DABS(PN) .GT. 1.D-10 ) THEN
            M = MAX0(m,MAX(II))
            DO 25 J = 1,M
 25           P(J,I) =P(J,I) - PN*P(J,II)
          END IF
	end if
6     CONTINUE
      pn = 1.d0/sqrt(quadr(i,i,0))
      if (p(4,i) .lt. 0.d0) pn = -pn
      do 16 j = 1,m
	p(j,i) = pn*p(j,i)
16     continue
      print *, ' expansion for ', n(i),l(i)
      print '(6f12.8)' ,(p(j,i),j=1,ns)
9     CONTINUE
*
*     .. improve estimates obtained from screened hydrogenics
*
      call improve(ind)
*
      WRITE(3,14)
14    FORMAT(/// 8X,18HINITIAL ESTIMATES  //10X,2HNL,
     1   4X,5HSIGMA,6X,5HE(NL),4X,9HFUNCTIONS//)
*
      DO 15 I = 1,NWF
      K = IND(I) + 2
      IF ( IND(I) .EQ. -2 ) THEN
           TITLE = ' SCALED '//ATM(I)//TRM(I)
        ELSE IF (IND(I) .EQ. 0) THEN
           TITLE = ' SCREENED HYDROGENIC'
        ELSE
           TITLE = ' UNCHANGED'
      END IF
17    WRITE(3,19) EL(I),S(I),E(I,I),TITLE
19    FORMAT(9X,A3,F9.2,F11.3,3X,A24)
15    CONTINUE
      ec = 0.d0
      IF (iuf .ne. 0) close(unit=iuf)
*      print *, ' Check orthogonlity'
*      a11 = quadr(1,1,0)
*      a12 = quadr(1,2,0)
*      a22 = quadr(2,2,0)
*      print *, a11, a12, a22
*      print *, '1s', '2s'
*      print '(6f12.8)', (p(ii,1),ii=1,ns)
*      print '(6f12.8)', (p(ii,2),ii=1,ns)
      END
*     ------------------------------------------------------------------
*               x c h
*     ------------------------------------------------------------------
*
*       this routine computes functions associated with the exchange
*   function for the i'th radial equation,  including   contributions
*   from   the  interactions. 
*
      SUBROUTINE xch(i,rel,xv)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER (nt=100,kx=8,ktx=3*kx-2,nwd=10)
      DOUBLE PRECISION xv(nt)
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rkm(nt,nt,kx,kx),lh,krk

      COMMON /SLATER/fc1(nt,kx),fc2(nt,kx),yk(nt),fyk(nt,kx),
     :        dyk(ktx,nt),ipvtd(nt),fc(nt,kx),ic1,ic2,iy1,iy2,ky
*
      LOGICAL kfirst, REL
      DOUBLE PRECISION xr(nt,kx)
*---------------------------------------------------------------------
*       
        nint = ns-ks+1
*       .. clear xr
	do 15 m = 1,ks
	  do 15 ii = 1,ns
	    xr(ii,m) = 0.d0
   15   continue
        if (ic1 .ne. i) then
           call yval(0,0,nt,kx,ks,nint,gr,bsp,p(1,i),fc1)
	   ic1 = i
        end if
*
*       .. compute contribution from each k
	do 20  k = 0,2*lmax
	  ykn   = 0.0d0
	  kfirst = .true.
	  do 30 j = 1,nwf
	    if (j .ne. i) then
	      if (iabs(l(i)-l(j)) .le. k .and. k.le. l(i)+l(j)) then
		if (mod(k,2) .eq. mod(l(i)+l(j),2)) then
                  cc = 2.d0*b(i,j,k)
		  print *, 'i,j,k,b(i,j,k)',i,j,k,cc
	          if (abs(cc) .gt. 1.d-10) then
		    if ( kfirst ) then
		      kfirst = .false.
*                     .. clear the array yk
	              do 45 ii=1,ns
	                yk(ii)=0.d0
 45 	              continue
*
*                     set up and factor the differential operator
*
                      call facdyk(nt,kx,ktx,ks,k,nint,t,bsp,bspd,
     :                            db2,rm2,ipvtd,dyk)
*   
		    end if
	            if (ic2 .ne. j) then
	              call yval(0,0,nt,kx,ks,nint,gr,bsp,p(1,j),fc2)
	              ic2 = j
	            end if
		    do 50 m = 1,ks
		      do 50 ii = 1,nint
		        fyk(ii,m) = cc*fc1(ii,m)*fc2(ii,m)
 50                 continue
	          end if
*
*                 add a relativistic correction if rel=TRUE
*
	          if (rel) then
	            c = (2*k+1)*fine
                    do 60 m = 1,ks
	              do 60 ii = 1,nint
	                xr(ii,m)= xr(ii,m)+ c*fyk(ii,m)
  60                continue
	          end if
*
* set up the array xk (spline approx. of SUM ()(1/r)fc1(r)fc2(r)
*
  	          do 70 m=1,ks
	            do 70 ith=1,ks
	              do 70 ii = 1,nint
         	        irow=ii+ith-1
	                yk(irow)=yk(irow)+
     :                      grw(ii,m)*grm(ii,m)*bsp(ii,m,ith)*fyk(ii,m)
 70	          continue
*
  	          const=-(2*k+1.d0)
	          do 80 ii=1,ns
	            yk(ii)=const*yk(ii)
 80	          continue
	          yk(1)=0.d0
	          yk(ns) = ykn
*
* solve the matrix equation
*
		  print *, 'yk before dgbsl'
		  print '(6f12.8)', (yk(ii),ii=1,ns)
	          call dgbsl(dyk,ktx,ns,ks-1,ks-1,ipvtd,yk,0)
                  call yval(0,0,nt,kx,ks,nint,grm,bsp,yk,fyk)
		  do 90 m = 1,ks
		    do 90 ii = 1,nint
		      xr(ii,m) = xr(ii,m) + fyk(ii,m)*fc2(ii,m)
 90               continue
		end if
	      end if
	    end if
 30       continue
*
 20    continue
*
*  evaluates the function x(r) at all the gaussian points
*
	do 100 m = 1,ks
	  do 100 ii = 1,nint
	    xr(ii,m) = xr(ii,m)*grw(ii,m)*grm(ii,m)
 100    continue
	call vinty(nt,kx,ks,nint,bsp,xr,xv)
*
	end
