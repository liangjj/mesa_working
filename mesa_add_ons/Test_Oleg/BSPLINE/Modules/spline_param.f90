!====================================================================
    MODULE spline_param
!====================================================================
!
!   contains basic spline parameters
!
!--------------------------------------------------------------------

    IMPLICIT NONE
    SAVE

    INTEGER(4) :: grid_type = 0   !  type of grid
    INTEGER(4) :: nug = 99 
    Character(40) :: AF_grid = 'knot.dat'

    INTEGER(4) :: ks = 8  !   order of B-splines
    INTEGER(4) :: ns = 0  !   number of splines
    INTEGER(4) :: nv = 0  !   number of intervals ( = ns-ks+1 )
    INTEGER(4) :: ml = 0  !   number of intervals from 0 to 1 (=1/h)
    INTEGER(4) :: me = 0  !   number of intervals in the exponential region

    REAL(8) :: h = 0.25    !   initial step in the knot sequence for z*r
    REAL(8) :: hmax = 0.25 !   maximum step, t(ns+1) - t(ns) 
    REAL(8) :: rmax = 1.00 !   border radius, t(ns+1)

    END MODULE spline_param
