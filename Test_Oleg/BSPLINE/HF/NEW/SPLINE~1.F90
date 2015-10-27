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
