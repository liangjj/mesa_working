!=========================================================================
    SUBROUTINE define_grid (z)
!=========================================================================
!
!   gets input data for the grid and sets up the knots for spline
!
!   SUBROUTINE contained:
!       getinput
!       mkgrid
!       print_grid
!
! -------------------------------------------------------------------------

    USE spline_param
    USE spline_grid

    IMPLICIT NONE
    REAL(8), INTENT(OUT) :: z

    ! .. Local variables

    INTEGER ::  nt

    ! .. get input data for the grid

    CALL getinput

    ! .. set up the knots for spline

     CALL mkgrid

    ! .. print grid in the file 'knot.dat'

    Call print_grid

    CONTAINS


!=======================================================================
      SUBROUTINE getinput
!=======================================================================
!  .. gets input data for the grid

      IMPLICIT NONE

      INTEGER :: nu
      Logical :: EX

      Inquire(file='knot.dat',exist=EX)

      IF(EX) THEN           ! read the grid parameters from file

       nu=99
       Open(nu,file='knot.dat')
       read(nu,*) ks
       read(nu,*) ns
       read(nu,*) z
       read(nu,*) h
       read(nu,*) hmax
       read(nu,*) rmax
       close(nu)

      ELSE                  ! read the grid parameters from terminal

       PRINT *, 'the following parameters are required for input:'
       PRINT *, 'real::    z for nuclei charge'
       PRINT *, 'real::    h for space step-size, of form 2**(-n) '
       PRINT *, 'real::    hmax for the largest space step allowed '
       PRINT *, 'real::    rmax for the maximun r of grid '
       PRINT *, 'integer:: ks for the order of B-spline  '
       PRINT *
       PRINT *, 'Enter z, h, hmax, rmax, ks (all on one line) '
       READ  *,  z, h, hmax, rmax, ks
       PRINT *
     END IF

    END SUBROUTINE getinput




!=====================================================================
    SUBROUTINE mkgrid
!=====================================================================

!   sets up the knots for spline

      IMPLICIT NONE

      ! .. Local variables

!      INTEGER, INTRINSIC:: NINT
      INTEGER:: n, i, m, me1, me2
!      REAL(KIND=8), INTRINSIC:: LOG, MAX
      REAL(KIND=8):: hp1, h2, tmax, tx, h0

      ! .. determine ml, the number of distinct points from 0 to 1

       h0 = h
       ml = 1.d0/h + 0.1
       h = 1.d0/ml
       hp1 = 1.d0 + h

      ! .. determine tmax

      tmax = z*rmax
      hmax = z*hmax

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

      ALLOCATE (t(nt)); t = 0.d0

      DO i = ks+1, ks+ml
        t(i) = t(i-1) + h0
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
      hmax = hmax/z

    END SUBROUTINE mkgrid


!=======================================================================
      SUBROUTINE print_grid
!=======================================================================

!  .. print grid parameters in 'knot.dat'

      IMPLICIT NONE

      INTEGER :: nu,i

      nu=99
      Open(nu,file='knot.dat')
       rewind(nu)
       write(nu,'(i12,a)') ks,' ==>  order  of splines (ks)'
       write(nu,'(i12,a)') ns,' ==>  number of splines (ns)'
       write(nu,'(f12.5,a)') z,' ==>  nuclear charge (z)'
       write(nu,'(f12.5,a)') h,' ==>  step size from 0 to 1 (h for z*r, = 1/2^n)'
       write(nu,'(f12.5,a)') hmax,' ==>  maximum step size (hmax for r)'
       write(nu,'(f12.5,a)') rmax,' ==>  maximum r (rmax)'
       write(nu,'(a)') '***'
       write(nu,'(5f12.5)') (t(i),i=1,ns+ks)
       write(nu,'(a)') '***'
       write(nu,'(i12,a)') ml,' ==>  number of distinct knots from 0 to 1 (=1/h)'
       write(nu,'(i12,a)') me,' ==>  number of knots in the exponential region '
       write(nu,'(f12.5,a)') (1/hmax)**2,' ==>  max k^2 (Ry) '

      close(nu)

    END SUBROUTINE print_grid

  END SUBROUTINE define_grid

