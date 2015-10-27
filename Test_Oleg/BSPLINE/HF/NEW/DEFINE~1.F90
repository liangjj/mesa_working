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
