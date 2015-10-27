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

