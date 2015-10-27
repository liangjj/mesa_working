!======================================================================
  PROGRAM Spline_HF
!======================================================================
!
!  This program computes the radial functions for simple Hartree-Fock
!  cases.
!
!   SUBROUTINE called:
!       get_case
!       define_grid
!       define_spline
!       define_atomic
!       obtain_estimates
!       optimize_orbitals
!       properties
!
!   Date: 01/21/1999
!                                                                
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    REAL(KIND=8) :: z

    ! .. get data about the problem to be solved
    CALL get_case(z)

    ! .. sets up positions of grid points: t
    CALL define_grid (z)

    ! .. initializes the values of the spline and its derivatives
    ! .. and evaluates the spline arrays (operators in spline basis)
    ! .. which are defined in the MODULES spline_grid and spline_galerkin
    CALL define_spline(z)

    ! ... initialize that atomic physics environment
    CALL define_atomic

    ! .. obtain initial estimates for orbitals
    CALL obtain_estimates

    ! .. optimize the orbitals for a stationary solution
    CALL optimize_orbitals

    ! .. compute atomic properties
    CALL properties

  END PROGRAM slater
!=======================================================================
  MODULE atomic_state 
!=======================================================================
!   This module defines the parameters for the problem to be solved
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    SAVE

    ! atom variables
    REAL(KIND=8) ::z
    INTEGER :: nclosd, nwf, nit
    CHARACTER(LEN=2) :: atom
    CHARACTER(LEN=3) :: term
    CHARACTER(LEN=32) :: config

    ! orbital variables
    CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE :: el
    INTEGER, DIMENSION(:), ALLOCATABLE :: n, l, max, meth, ind
    INTEGER :: lmax
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: p
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: az, dpm, s, sum

    ! energy expression variables
    INTEGER :: kmax
    INTEGER, DIMENSION(:), ALLOCATABLE :: ijptr
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: coef
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: e

    CONTAINS

    !===================================================================
      SUBROUTINE allocate_atomic_state
    !===================================================================
    !   This program allocates arrays associated with number of orbitals
    !-------------------------------------------------------------------
        IMPLICIT NONE
        
	ALLOCATE( n(nwf), l(nwf), max(nwf), ind(nwf), s(nwf), meth(nwf))
	ALLOCATE( az(nwf), dpm(nwf), s(nwf), sum(nwf) )
        ALLOCATE( e(nwf,nwf) )

      END SUBROUTINE allocate_atomic_state

    !===================================================================
      SUBROUTINE allocate_orbital_array
    !===================================================================
    !   This program allocates arrays associated with number of orbitals
    !   and the spline expansion
    !-------------------------------------------------------------------
        USE spline_param
        IMPLICIT NONE
        
	ALLOCATE( p(ns,nwf) )
        ALLOCATE( ijptr(nwf-nclosd, nwf-nclosd) )
     
      END SUBROUTINE allocate_orbital_array
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
 
  CONTAINS
  !=====================================================================
    SUBROUTINE bwint(lc,lo)
  !=====================================================================
  !
  !   Compute the spin-orbit parameter using the Blume and Watson method
  !
  !----------------------------------------------------------------------
  !
  
    IMPLICIT NONE
      COMMON/blume/coefn2(4),coefnk(4),coefvk(4)
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
  
  END SUBROUTINE bwzeta
!=========================================================================
  SUBROUTINE define_grid (z)
!=========================================================================
!  
!   gets input data for the grid and sets up the knots for spline 
!
!   SUBROUTINE contained:
!       getinput
!       mkgrid
!
!   calling sequence:
!       define_grid
!       -----------
!         //    \\
!     getinput mkgrid
!
! -------------------------------------------------------------------------   
!
    USE spline_param
    USE spline_grid

    IMPLICIT NONE
    REAL(KIND=8), INTENT(OUT) :: z

    ! .. Local variables
    INTEGER ::  nt
    REAL(KIND=8)::  hmax,rmax           
  
    ! .. get input data for the grid
    CALL getinput

    ! .. set up the knots for spline
    CALL mkgrid

    CONTAINS

    !=====================================================================
    SUBROUTINE getinput
    !=====================================================================
    !   gets input data for the grid
    !---------------------------------------------------------------------   
    !
      IMPLICIT NONE
      PRINT *, 'the following parameters are required for input:'
      PRINT *, 'real::    h for space step-size, of form 2**(-n) '
      PRINT *, 'real::    hmax for the largest space step allowed '
      PRINT *, 'real::    rmax for the maximun r of grid '
      PRINT *, 'integer:: ks for the order of B-spline  '
      PRINT *
      PRINT *, 'Enter  h, hmax, rmax, ks (all on one line) '
      READ *, h, hmax, rmax, ks
      PRINT *   
    END SUBROUTINE getinput

    !=====================================================================
    SUBROUTINE mkgrid
    !=====================================================================
    !   sets up the knots for spline
    !---------------------------------------------------------------------   
    !
      IMPLICIT NONE
      ! .. Local variables
      ! .. INTEGER:: ml, me
      INTEGER, INTRINSIC:: NINT
      INTEGER:: n, i, m, me1, me2
      REAL(KIND=8), INTRINSIC:: LOG, MAX
      REAL(KIND=8):: hp1, h2, tmax, tx

      ! .. determine ml, the number of distinct points from 0 to 1
      ml = NINT(1.d0/h)
      h = 1.0d0/ml
      hp1 = 1.d0 + h
 
      ! .. determine tmax
      tmax = z*rmax
 
      ! .. determine final point of "exponential" grid
      ! .. me: number of points from 1 to (1+h)**me
      ! .. m:  number of points from (1+h)**me to tmax
      me1 = MAX(0.0d0, LOG(hmax/h)/LOG(hp1)+1.d0)
      me2 = LOG(tmax)/LOG(hp1)+1

      IF ( me2 <= me1 ) THEN
        me = me2
        m = 0
      ELSE 
        me = me1
        tx = hp1**me
        h2 = h*tx/hp1
        m = NINT((tmax-tx)/h2)
      END IF
      n = ml + me + m + ks -1
      ns = n
      nv = ns - (ks -1)
      nt = ns + ks
      
      ! .. establish the grid for z*r     

      ALLOCATE (t(nt))

      t(1:ks) = 0.d0

      DO i = ks+1, ks+ml
        t(i) = t(i-1) + h
      END DO

      DO i = ks+ml+1, ks+me+ml
        t(i) = t(i-1)*hp1
      END DO

      DO i = ks+me+ml+1, n+1
        t(i) = t(i-1) + h2
      END DO
      t(n+2:nt) = t(n+1)

      ! .. scale the values to the R variable
      t = t/z
    END SUBROUTINE mkgrid
   
  END SUBROUTINE define_grid
!======================================================================
  SUBROUTINE define_spline(z,l)
!======================================================================
!                                                                  
!   initializes the values of the spline and its derivatives
!   and evaluates the spline basic arrays (elementary operators in 
!   spline basis).
!	                                                         
!   SUBROUTINE called:
!       gauss 
!       allocate_memory
!       initvb
!       initas
!       hlm
!                                                                
!   calling sequence:
!                          define_spline   
!                   ------------------------- 
!                  / |      |       |      | 
!                 /  |   initvb  initas    |
!                /   |      |     /  \     |  
!           gauss    |   vbsplvd mdb mrm  hlm  
!                    |      ||                
!       allocate_memory  vbsplvb         
!         
!----------------------------------------------------------------------
!
    USE spline_param

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(ks) :: gx, gw
    REAL(KIND=8) :: z

    ! .. initializes variables for gaussian integration
    CALL gauss(ks,gx,gw)

    ! .. allocates space for arrays defined in MODULE spline_grid
    ! .. and spline_galerkin.
    CALL allocate_memory

    ! .. initializes the values of the spline and its derivatives
    CALL initvb(gx,gw)

    ! .. initializes the spline array (operators in spline basis)
    CALL initas

  END SUBROUTINE define_spline

 
!=======================================================================
  SUBROUTINE get_case(z) 
!=======================================================================
!  Get information from the user about the problem to be solved
!----------------------------------------------------------------------
!
    USE atomic_state
    USE periodic_table
    IMPLICIT NONE

    CHARACTER(LEN=1) :: ans
    LOGICAL :: ok

    write(*,'(A)',ADVANCE='NO') 'Atomic Number                         : '
    read(*,*) z
    atom = atom_label(nint(z))
    
    DO
      ! obtain closed shells
      call closed_shells
  
      ! obtain configuration
      write(*,'(A)',ADVANCE='NO') 'Configuration [ Ex. 2s(1)2p(3) ]      : '
      read(*,*) config
  
      ! obtain term
      write(*,'(A)',ADVANCE='NO') 'Term [ Ex. 3P or AV ]                 : '
      read(*,*) term
  
      ! determine the list of electrons
      call find_electrons
      IF (ok) EXIT
    END DO
    
    ! determine orbitals to be varied
    DO
      write(*,'(A,I3,A)') 'There are',nwf,' orbitals as follows:'
      write(*,'(4x,18(1X,A3))') el
      write(*,'(A/A)',ADVANCE='NO') 'Indicate those to be varied', &
                                '[Ex. all/none/=i/comma delimited list/h]:'
      call orbitals_varied
      if (ok) EXIT
    END DO

    ! determine orbital parameters
    write(*,'(A)',ADVANCE='NO') 'Default orbital parameters  (y/n)    : '
    read(*,*) ans
    call orbital_parameters(nwf,el,nclosd,ans)
    
    ! get the energy expression
    call get_energy_expression
   
    CONTAINS 

    !===================================================================
      SUBROUTINE closed_shells
    !===================================================================
    !  Get information from the user about closed shells
    !------------------------------------------------------------------
    !
      IMPLICIT NONE
      INTEGER :: i,ii, core_type

      i = 0
      DO
	IF (number_of_electrons(i+1) > z) then
	  EXIT
        ELSE
          i = i + 1
        END IF
      END DO

      Write(*,*) 'Possible cores are:
      Write(*,'(3X,I3,A,2X,A)') (ii,'.', shells(ii),ii=1,i)
      Write(*,*,ADVANCE='no') 'Enter the number for your selection: '
      Read(*,*) core_type
      i = 0
      DO
	IF (shell_number(i+1) > number_or_electrons(core_type)) then
	  EXIT
        ELSE
          i = i + 1
        END IF
      END DO
      nclosed = i
     
    END SUBROUTINE closed_shells

    !===================================================================
      SUBROUTINE find_electrons
    !===================================================================
    !  Analyze the configuration and report the list of electrons
    !------------------------------------------------------------------
    !
      IMPLICIT NONE
      INTEGER :: i,is,ipl,ipr,iel
      CHARACTER(LEN=3) :: ell
      CHARACTER(LEN=80) :: sconfig

      ok = .true.
      write(*,'(A,F6.0,A,A)')'Atom         : ', z,atom
      write(*,'(A,18(1X,A2)) 'Closed Shells: ',(shell_order(1:nclosd))
      Write(*,*)             'Configuration: ', trim(config), '  ',term
      ! count the number of shells in the configuration
      i = 0
      is = 1
      DO
       ipl = index(config(is:),'(')
       IF (ipl > 0) then 
         ! there was a shell
         IF (ipl-is > 3) then
           write(*,*) 'ERROR: Electron label ', config(is:ipl-1), &
                      ' has more than 3 characters'
           ok = .false.
           RETURN
         END IF
	 ell = config(is:ipl-1)
	 ! check this is not in closed shells 
         DO iel = 1,nclosd
           IF (ell == shell_order(iel)) then
             write(*,*) 'ERROR: Electron', ell, ' appears both in ', &
                        config, 'and the closed shells'
             ok = .false.
             RETURN
           END IF
         END DO
	 ipr= index(config(ipl+1:),')')
         ! check for matching parentheses 
         IF (ipr = 0 ) then
           write(*,*) 'ERROR: Matching parentheses not found in configuration'
           ok = .false.
           RETURN
         ELSE IF (ipr-ipl .gt. 5) then
           write(*,*) 'ERROR:'                                       &
             'Occupation number should not exceed 4 characters'
           ok = .false.
           RETURN 
         END IF
	 i = i+1
	 is = ipr+1
       ELSE
	 EXIT
       END IF
      END DO
      nwf = nclosd+i
 
      ! at this point the orbital data should be correct
      call allocate_atomic_state

      ! set the list of orbitals
      el(1:nclosd) = shell_order(1:nclosd)
      el(1:nclosd) = ADJUSTR(el(1:nclosd))
      is = 1
      call reform(TRIM(ADJUSTL(config)), sconfig)
      read(sconfig,'(8(1x,a3,6x)))') el(nclosd+1:nwf)
 
    END SUBROUTINE find_electrons
      
  !===================================================================
    SUBROUTINE orbitals_varied
  !===================================================================
  !  Determine which orbitals are to be varied
  !------------------------------------------------------------------
  !
     IMPLICIT NONE
     CHARACTER(LEN=72) :: string
     
     ok = .true. 
     read(5, '(a)') string
     if (string(1:1) .eq. 'h' .or. string(1:1) .eq. 'H') then
       call help(1)
       ok = .false.
       RETURN
     else if (string(1:3) .eq. 'all' .or. string(1:3) .eq. 'ALL') then
       nit = nwf
     else if (string(1:4).eq.'none' .or. string(1:4).eq.'NONE') then
       nit = 0
     else if (index(string,'=') .ne. 0) then
       j = index(string,'=')
       jj = index(string,' ')
       if ( jj .eq. j+2) then
	 read(string(j+1:j+1),'(i1)') nit
       else if (jj .eq. j+3) then
	 read(string(j+1:j+2),'(i2)') nit
       else
	 write(iscw,'(a,a)') ' nit must be specified by one or two' &
              'digits immediately following = sign (no blanks)
	 ok = .false.
         RETURN
       end if
     else
       nit = 0
       j = 1
       DO
         next = index(string(j:),',')
         ! search for last electron LABEL which need not be followed by a comma 
         if (next .eq. 0 .and. string(j:j+2) .ne. '   ') &
             next = index(string(j+1:),' ') + 1
         if (next .ge. 1) then
            if (next .eq. 4) then
               el1 = string(j:j+2)
            else if (next .eq. 3) then
               el1 = ' '//string(j:j+1)
            else
               write(iscw,*) &
                 ' ERROR: electron LABELs must be separated by commas;'
               write(iscw,*)' each LABEL must contain 2 or 3 CHARACTERs'
               ok = .false.
               RETURN
            end if
            call reord(el,el1,nwf,ierr)
            if (ierr .eq. 0) then
               ! no error has occurred
               nit = nit + 1
               j = j + next
               if (next > 72) EXIT
            else
               write(*,*) 'ERROR: Electron ', el1, ' not found in the list'
               write(*,*) 'Case and position of imbedded blanks must match'
               ok = .false.
               RETURN
            end if
         END DO
         end if
      end if
    END SUBROUTINE orbitals_varied
      
  !===================================================================
    SUBROUTINE orbital_parameters
  !===================================================================
  !  Set orbital parameters
  !------------------------------------------------------------------
  !
     IMPLICIT NONE
     INTEGER :: i
     REAL(KIND=8) :: ss
      
     ! Set occupation number of non-core orbitals
     Read(sconfig,'(8(5x,F4.0,1X)') (sum(nclosd+1:nwf)
     ! Determine nl values and related quantities
     
     ss = 0
     lmax = 0
     DO i = 1,nwf
       l(i) = lval(el(i)(3:3)
       read(el(i)(1:2), '(I2)')  n(i)
       IF (i <= nclosd) sum(i) = 2*( 2*l(i) + 1)
       s(i) = ss + sum(i)/2
       ss = ss + sum(i)
       meth(i) = 1
       acc(i) = 0
       lmax = max(lmax,l(i))
     END DO

     ! For HF calculations we set
     strong = .false.

     kmax = 2*lmax

    END SUBROUTINE orbital_parameters
     
  !===================================================================
    SUBROUTINE get_energy_expression
  !===================================================================
  !  Obtain the data that defines the HF equations from an
  !  energy expression stored in a coef list and an ijptr array.
  !------------------------------------------------------------------
  !
     IMPLICIT NONE
     INTEGER :: i,j, ip
     REAL(KIND=8) ::C, CA, CB`

     ! Find the size of the coef array
     ip = 0
     DO i=1,nwf
       DO j= 1,nwf
	 ! direct contribution
         ip = ip +1
         if (i == j ) ip = ip + l(i)
         ! exchange contribution
         ip = ip + min(l(i),l(j))
       END DO
     END DO

     ALLOCATE ( coef(ip) )

     ip = 0
     do i = nclosd+1,nwf
       isumi = sum(i)
       dsumi = sum(i) - isumi
       do j = nclosd+1,nwf
         isumj = sum(j)
         dsumj = sum(j) - isumj
         if ( i .ne. j) then
           c = sum(j)
           if (dsumi .ne. 0.d0 .and. dsumj .ne. 0.d0) &
             c = (dsumi*(isumi+1)*isumj +a            &
                  dsumj*(isumj+1)*isumi)/sum(i)
           else
             c = sum(i) - 1.d0
             if (dsumi .ne. 0.d0)                     &
               c = (isumi*(sum(i)+dsumi-1))/sum(i)
           end if

           ijptr(i-nclosd,j-nclosd) = ip

!           ...direct contribution

           do k = 0,2*min0(l(i),l(j)),2
             ip = ip + 1
             if (k .eq. 0) then
               coef(ip) = c
             else if (i .eq. j) then
               coef(ip) = -c*ca(l(i),k)
             end if
           end dO

!           ... exchange contribution

           if (i .ne. j) then
             do k = abs(l(i)-l(j)),l(i)+l(j),2
               ip = ip + 1
               coef(ip) = -c*cb(l(i),l(j),k)
             END Do
           end if
         END Do
       END Do

       CALL enexpR

    END SUBROUTINE get_energy_expressioN
     
  !===================================================================
    SUBROUTINE enexpr
  !===================================================================
  ! Determine the deviations from the average energy for the term
  ! dependent energy expression for the following cases:
  !      i) an open p- or d-shell
  !     ii) a single electron or hole, any l
  !    iii) an s-electron and a single electron, any l
  !     iv) an s-electron and an open p- or d-shell
  !      v) an open p-shell and a single electron, any l
  !------------------------------------------------------------------
  !
      IMPLICIT NONE

      iNTEGER, DIMENSION(5) ::  sumtab(5)
      INTEGER, DIMENSION(11) :: partab,ptrtab,plval
      INTEGER, DIMENSION(54) :: ltaB
      INTEGER, DIMENSION(2)  :: nos(2)
      INTEGER :: pacval,sp,ps1,ps2
 
      ... fint, gint1, and gint2 are coefficients of polynomials
          in l, tabulated by slater,
 
 
      cHARACTER(LEN=1) :: parch(11)

      !... coefficients of f2 integrals for p(n)l(1) configurations
      INTEGER, DIMENSION(3,54) :: fint = &
        (/2,-1,0,-4,-4,3,2,5,3,2,-1,0,-4,-4,3,2,5,3,             &
         -2,1,0,4,4,-3,-2,-5,-3,-2,1,0,4,4,-3,-2,-5,-3,          &
         4,-2,0,-2,-11,6,-4,-4,15,-2,7,15,4,10,6,0,0,0,          &
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  &
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  &
         2,-1,0,-4,-4,3,2,5,3,2,-1,0,-4,-4,3,2,5,3,              &
         -4,2,0,2,11,-6,4,4,-15,2,-7,-15,-4,-10,-6,0,0,0,        &
         -2,1,0,4,4,-3,-2,-5,-3,-2,1,0,4,4,-3,-2,-5,-3/)
 
      !... coefficients of g(l-1) integrals
      INTEGER, DIMENSION(3,54) :: gint1 `= &
       (/-10,5,0,2,11,-6,2,-1,-6,14,-7,0,2,-13,6,2,-1,6,               &
         -8,4,0,-8,4,0,4,10,0,10,-5,0,10,-5,0,4,-8,0,                  &
         -8,4,0,-2,13,-6,2,11,-12,4,4,-12,4,-2,-12,0,0,0,-6,3,0,10,-5,0,&
         -6,3,0,-6,3,0,-4,8,3,0,12,3,6,9,0,-6,3,0,6,21,12,12,12,-27,    &
         12,-6,27,6,-15,-12,-6,3,0,0,6,-3,0,0,-3,6,-3,0,0,-6,3,12,6,3,  &
         -4,2,0,-4,2,0,-4,2,0,-4,2,0,14,11,-9,14,-7,-9,-4,-2,0,-4,2,0,  &
         -2,7,3,2,11,3,8,8,0,0,0,0,-2,1,0,-2,1,0,-2,1,0,-2,1,0,-2,1,0,  &
         22,13,0/)
 
      !... coefficients of g(l+1) integrals
      INTEGER, DIMENSION(3,54) :: gint2 `= &
       (/2,5,-3,2,-7,-15,-10,-25,-15,2,5,9,2,17,21,14,35,21,       &
         4,-2,-6,-8,-20,-12,-8,-20,-12,4,16,12,10,25,15,10,25,15,  &
         4,10,0,4,4,-12,2,-7,-21,-2,-17,-21,-8,-20,-12,0,0,0,      &
         -6,-15,-9,10,25,15,6,3,-3,0,-12,-9,-4,-16,-9,-6,-15,-9,   &
         -6,-15,-9,6,27,9,12,30,-9,12,12,-27,6,-9,-27,-6,-15,-9,   &
         0,0,-3,0,-6,-9,-6,-15,-9,12,18,9,0,6,9,6,15,9,            &
         -4,-10,-6,-4,-10,-6,-4,-10,-6,14,35,12,14,17,-6,-4,-10,-6,&
         8,8,0,2,-7,-6,-2,-11,-6,-4,-10,-6,-4,-10,-6,0,0,0,        &
         -2,-5,-3,-2,-5,-3,-2,-5,-3,22,31,9,-2,-5,-3,-2,-5,-3/)
 
      !.. encoded term value --   s = ltab/10
      !                       lterm = l + (ltab mod 10 - 5)
      !   example: ltab = 36 with l = 2  is 3f
 
      data ltab/36,35,34,16,15,14,46,45,44,26,25,24,27,26,25,24,
     :   23,25,55,35,37,36,35,34,33,17,16,15,14,13,36,35,34,16,15,14,
     :   46,45,44,26,25,24,27,26,25,24,23,25,36,35,34,16,15,14/
      data sumtab/1,4,7,10,11/
      data partab/2,3,1,1,4,2,2,3,1,1,2/
      data ptrtab/6,12,17,18,20,30,36,42,47,48,54/
      data plval/1,1,2,0,0,2,1,1,2,0,1/
      data parch/'p','p','d','s','s','d','p','p','d','s','p'/
 
       
      sl = term
      senor = ' '
 
      !  convert lowercase l symbol to uppercase
 
      if (sl(2:2).gt.'a' .and. sl(2:2).lt.'z')                      &
          sl(2:2) = char(ichar(sl(2:2)) + ichar('a') - ichar('a'))
 
      !  determine if fk or gk data needs to be input
 
      il = 0
      is = 0
      j = 1
      do I = nclosd+1, nwf
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
      END Do
      if (sl .ne. 'av' .and. sl .ne. 'av' .and. is+il.ne.0) then
         done = .false.
         c = 0.d0
         if (is+il .le. 2 .and. il .le. 1) then
            if (is .eq. 0 .and. il .eq. 1) then
	      Do
               call looktm(l(iil),sl,senor,sum(iil),ip,nsl)
	       if (nsl > 1) then
		  write(iscw,*)' ambiguous term: enter seniority'
		  read (5,'(a1)') senor
	       else
		  EXIt
	       end if
	      END Do
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
		 Do
                  write(iscw,'(a,a3,a,a3)')                          &
                       ' ambiguous l**n term: enter',slm, ' or ',slp 
	          read(5,'(a2)') sl
		  if (sl .eq. slm) then
		     ip = ipm
		     EXIt
		  else if (sl .eq. slp) then
		     ip = ipp
		     EXIt
		  else
	             write(iscw,*) ' term not allowed: re-enter'
		  end if
		 END Do
	       else 
		 Do
                  write(iscw,'(a,a)') ' ambiguous l**n parent term:', &
                         'enter term and seniority'
	          read(5,'(a2,a1)') sl, senor
	          call looktm(l(iil),sl,senor,sum(iil),ip,nsl)
	          if (nsl <> 1) then
		    write(iscw,'(a,a3,a,a3,a)') ' allowed terms are ',  &
                          slm, ' or ', slp,' plus seniority'
		  ELSE
		     EXIt
	          end if
		 END Do
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
            if (((l(nos(1)).eq.1).and.(sum(nos(2)).eq.1.d0)).or.    &
                ((l(nos(2)).eq.1).and.(sum(nos(1)).eq.1.d0))) then
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
		 Do
                  write(iscw,*)' enter the sl term for p(n) subshell'
                  read (5,'(a)')psl
                  if (psl(2:2).gt.'a'.and.psl(2:2).lt.'z')           &
                         psl(2:2)=char(ichar(psl(2:2))+ichar('a')    &
                                                      -ichar('a'))  
                  ps1=ichar(psl(1:1))-ichar('0')
                  ps2=lval(psl(2:2))
                  call lookup(plval,iptr1,iptr2,ind,nomach,ps2)
                  if (.NOT.(nomach.ne.1).and.(partab(ind).ne.ps1)) EXIT
                 END Do
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
               val1=((fint(1,ind)*lv+fint(2,ind))*lv+fint(3,ind))    &
                        /(5.d0*(2*lv-1)*(2*lv+3))
               val2=((gint1(1,ind)*lv+gint1(2,ind))*lv+gint1(3,ind)) &
                        /(2.d0*(2*lv+1)*(2*lv-1)**2)
               val3=((gint2(1,ind)*lv+gint2(2,ind))*lv+gint2(3,ind)) &
                        /(2.d0*(2*lv+1)*(2*lv+3)**2)
 
               !  add contributions from between p-subshell and l-electron
 
               call add(val1,2,np,nl,.true.)
               call add(val1,2,nl,np,.true.)
               call add(val2,lv-1,np,nl,.false.)
               call add(val2,lv-1,nl,np,.false.)
               call add(val3,lv+1,np,nl,.false.)
               call add(val3,lv+1,nl,np,.false.)
 
               !... add deviations for p-subshell
 
               call looktm(1,psl,' ',sum(np),ip,nsl)
               call dev(np,1,sum(np),ip,done)
            else
               done = .false.
            end if
         end if
      else
         done = .true.
      end if
     END SUBROUTINE enexpr
 
  !=====================================================================
    SUBROUTINE looktm( l, sl, sen, q, ip, nsl)
  !=====================================================================
  !
  !    add the deviations to the average energy for a partially filled
  !      p- or d- shell
  ! 
  !--------------------------------------------------------------------
  !
    IMPLICIT NONE
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
         DO
         if (sl .eq. terms(i)(1:2)) then
            if (sen .eq. ' ' .or. sen .eq. terms(i)(3:3)) then
	       nsl = nsl + 1
	       ip = i
            end if
         end if
         i = i+1
         if (i > iend) EXIT
         END DO
      else if ( n .eq. 1 .and. sl(1:1) .eq. '2') then
	 nsl = 1
      end if
    END SUBROUTINE looktm
 
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
     END SUBROUTINE dev

   END SUBROUTINE get_case
!==================================================================
  SUBROUTINE initas
!==================================================================
!
!   Sets ( or Initializes ) the array in symmetric storage mode:
!
!       db1 --- matrix of integral <B_i,B'_j>
!       db2 --- matrix of integral <B_i,B"_j>
!       sb  --- matrix of integral <B_i,B_j>
!       r1  --- matrix of integral <B_i,r B_j>
!       rm1 --- matrix of integral <B_i,(1/r)B_j>
!       rm2 --- matrix of integral <B_i,(1/r^2)B_j>
!               where i=1,..,ns, j=1,...ks
!
!   SUBROUTINES called
!       mdb mrm
!
!   Calling sequence:
!       initas
!        /  \
!      mdb  mrm
!
!------------------------------------------------------------------
!
    USE spline_galerkin

    IMPLICIT NONE

    ! .. sets db1 --- matrix of integral <B_i,B'_j>
    CALL mdb(1, db1)  
 
    ! .. sets db2 --- matrix of integral <B_i,B"_j>
    CALL mdb(2, db2)  
 
    ! .. sets sb  --- matrix of integral <B_i,B_j>
    CALL mrm(0, sb)  

    ! .. sets r1  --- matrix of integral <B_i,r B_j>
    CALL mrm(1, r1)  

    ! .. sets rm1 --- matrix of integral <B_i,(1/r)B_j>
    CALL mrm(-1, rm1)  
 
    ! .. sets rm2 --- matrix of integral <B_i,(1/r^2)B_j>
    CALL mrm(-2, rm2)  

  END SUBROUTINE initas

!=======================================================================
  MODULE periodic_table 
!=======================================================================
!   This module lists the labels of the elements of the periodic table 
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    SAVE

    CHARACTER(LEN=2), DIMENSION(103) ::  atom_label = &
     (/' H', 'He', 'Li', 'Be', ' B', ' C', ' N', ' O', ' F', 'Ne', &
       'Na', 'Mg', 'Al', 'Si', ' P', ' S', 'Cl', 'Ar', ' K', 'Ca', &
       'Sc', 'Ti', ' V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
       'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', ' Y', 'Zr', &
       'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
       'Sb', 'Te', ' I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', &
       'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
       'Lu', 'Hf', 'Ta', ' W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
       'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
       'Pa', ' U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', &
       'Md', 'No', 'Lw', '  ', '  ', '  ', '  ', '  ', '  ', '  '/)

    CHARACTER(LEN=2), DIMENSION(17) :: shell_order = &
     (/'1s', '2s', '2p', '3s', '3p', '3d', '4s', '4p', '4d', '5s', &
       '5p', '4f', '5d'  '6s', '6p', '5f', '7s', '7p'/
 
    INTEGER, DIMENSION(17) :: shell_number = &
     (/   2,    6,  10,   12,   18,   28,   30,   36,   46,   48,  &
         54,   68,  78,   80,   86,  100,  102,  108/)

    CHARACTER(LEN=*), DIMENSION(13) :: shells = &
     (/'1s(2)', '.. 2s(2)2p(6)', '.. 3s(2)3p(6)',                  &
       '.. 3d(10)', ' .. 4s(2)4p(6)', '4d(10)', '.. 5s(2)5p(6)',   &
       '.. 4f(14)', ' '.. 5d(10)', '.. 6s(2)6p(6)' &
       '.. 5f(14)', ' .. 7s(2)7p(6)'/)

    INTEGER, DIMENSION(13) :: number_of_electrons = &
     (/  2,   10,   18,   28,   36,   46,   54,    ` &
        68,   78,   86,  100,  108/)

  END MODULE periodic_table



!==========================================================================
  MODULE spline_galerkin
!==========================================================================
!  The SPLINE_GALERKIN module contains common arrays used in the
!  application of splines and the Galerkin method.
! -------------------------------------------------------------------------   
    IMPLICIT NONE
    SAVE

    ! .. spaces for initializing spline values and arrays
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:):: r1,rm1,rm2
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:):: sb, hl, db1, db2
  
    CONTAINS

    !======================================================================
      SUBROUTINE allocate_galerkin
    !======================================================================
    !  This program allocates space for the arrays used in the
    !  application of splines and the Galerkin method.
    !----------------------------------------------------------------------   
        USE spline_param
        IMPLICIT NONE
        ALLOCATE( r1(ns,ks),rm1(ns,ks),rm2(ns,ks) )
        ALLOCATE( sb(ns,ks), hl(ns,ks), db1(ns,ks),db2(ns,ks) )
      END SUBROUTINE allocate_galerkin

  END MODULE spline_galerkin
!==========================================================================
  MODULE spline_grid
!==========================================================================
!  The spline_grid module defines the values of splines at the gaussian
!  points defined by the intervals of a grid.  Included in the module
!  is the gaussian data for performing integrations on the grid.
! -------------------------------------------------------------------------   
    IMPLICIT NONE
    SAVE

    ! .. arrays for defining grid
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE:: t

    ! .. arrays for initializing spline values 
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: bs
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE:: bsp     
    REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE:: bspd  

    ! .. arrays for initializing gaussian data
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE:: gr, grm, grw

    CONTAINS

    !======================================================================
      SUBROUTINE allocate_grid
    !======================================================================
    !  This program allocates space of the arrays for initializing spline
    !  values and gaussian data in MODULE spline_grid
    !----------------------------------------------------------------------   
        USE spline_param
        IMPLICIT NONE
        ALLOCATE( bs(ks,ns), bsp(nv+1,ks,ks), bspd(nv+1,ks,ks,2) )
        ALLOCATE( gr(nv,ks),grm(nv,ks),grw(nv,ks) )
      END SUBROUTINE allocate_grid

  END MODULE spline_grid

!==========================================================================
  MODULE spline_param
!==========================================================================
!  The SPLINE_PARAM module defines the spline parameters
! -------------------------------------------------------------------------
    IMPLICIT NONE
    SAVE

    ! .. Here are the commonly used parameters
    INTEGER :: ks, ns, nv,   ml,me
    REAL(KIND=8) :: h

  END MODULE spline_param
!=======================================================================
  SUBROUTINE wavefn 
!=======================================================================
!
!       This routine initializes radial functions.  The file 'wfn.inp"
!  is read.  The first radial function with the same electron label is
!  scaled to the current Z as an initial estimate.  Otherwise, an 
!  estimate is formed as a screened hydrogenic function witn ZZ=Z-s(i)
!  and a Hartree-orbital computed, where previously defined orbitals
!  define the potential and the new orbital is constrained to be orthogonal
!  to existing orbitals.
!
!----------------------------------------------------------------------
!
    IMPLICIT NONE
!
      CHARACTER*3 EL1
      CHARACTER*6 AT,TT,ATM(NWD),TRM(NWD)
      CHARACTER*24 TITLE
      LOGICAL LD
!
      nint=ns-ks+1
      call facsb(nt,kx,ks,ns,sb,bs)
!
!  ***** READ THE WAVEFUNCTIONS
!
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
!
!  *****  SCALE RESULTS IF DATA IS FOR AN ATOM WITH A DIFFERENT Z
!
         e(I,I) = C*C*ETI
         DO 11 J = 1,M
            P(J,I) = C*PT(J)
11       CONTINUE
!
!  *****  SET REMAINING VALUES IN THE RANGE = 0.
!
         IF ( M .EQ. ns ) GO TO 12
         M = M +1
         DO 13  J=M,ns
13       P(J,I) = 0.d0
12       IND(I) = -2
      ENDIF
      GO TO 2
!
!  *****  SET PARAMTERS FOR ELECTRONS AND INITIALIZE FUNCTIONS
!
5     continue
      DO 9 I = 1,NWF
      IF (IND(I)) 7,8,9
!
!  ***** WAVE FUNCTIONS NOT FOUND IN THE INPUT DATA, SET IND = 0
!
7     IF ( IND(I) .EQ. -2 ) GO TO 4
      IND(I) = 0
      WRITE(iscw,27) EL(I)
27    FORMAT(8X,'WAVE FUNCTIONS NOT FOUND FOR ',A3)
!
!  *****  DETERMINE ESTIMATES OF THE WAVE FUNCTIONS BY THE SCREENED
!  *****  HYDROGENIC APPROXIMATION
!
8     continue
      CALL BHWF(N(I),L(I),Z-s(i),nt,kx,ks,nint,gr,grw,bsp,bs,w21,w22,
     :          P(1,I))
      M = NS-1
30    IF ( DABS(P(M,I)) .LT. 1.D-15 ) then
        P(M,I) = 0.d0
        M = M-1
        GO TO 30
      END IF
!31    MAX(I) = M+1
 31    MAX(I) = ns
       print *, ' nl expansion', el(i)
       print '(6f12.8)', (p(ii,i),ii=1,ns)
      e(i,i) = ((Z-s(i))/n(i))**2
       print *, 'After orthogonalizationn'
!
!  *****  ORTHOGONALIZE TO INNER FUNCTIONS
!
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
!
!     .. improve estimates obtained from screened hydrogenics
!
      call improve(ind)
!
      WRITE(3,14)
14    FORMAT(/// 8X,18HINITIAL ESTIMATES  //10X,2HNL,
     1   4X,5HSIGMA,6X,5HE(NL),4X,9HFUNCTIONS//)
!
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
!      print *, ' Check orthogonlity'
!      a11 = quadr(1,1,0)
!      a12 = quadr(1,2,0)
!      a22 = quadr(2,2,0)
!      print *, a11, a12, a22
!      print *, '1s', '2s'
!      print '(6f12.8)', (p(ii,1),ii=1,ns)
!      print '(6f12.8)', (p(ii,2),ii=1,ns)
      END
