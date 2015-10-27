!***********************************************************************
! FEDVR_Module
!**begin prologue     FEDVR_Module
!**date written       090219   (yymmdd)
!**revision date               (yymmdd)
!**keywords           DVR, FEDVR
!**
!**author             schneider, b. i.(nsf)
!**source             DVR Library
!**purpose            Driver module to calculate the FEDVR functions and matrix elements
!***                  in a FEDVR basis
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      FEDVR_Module
!***********************************************************************
!***********************************************************************
                           MODULE FEDVR_Module
                           USE FEDVR_Global
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Lobatto_Functions.f
!***begin prologue     Lobatto_Functions
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculate piecewise dvr functions
!***description        DVR polynomials are computed for a number                   
!***                   of coordinate systems.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Lobatto_Functions

  SUBROUTINE Lobatto_Functions(grid)
  IMPLICIT NONE
  TYPE (coordinates)            :: grid
!
  IF (keyword == 'cartesian'.or.keyword == 'legendre'.or.keyword == 'theta'   &
                            .or.keyword == 'hermite'.or.keyword == 'laguerre') THEN
      ALLOCATE ( grid%reg_pt_wt(1:nreg), grid%reg_poly(1:nreg) )
      Call LaGrange_Polynomials(grid,grid%reg_poly)
  ELSE IF(keyword =='spherical') THEN
      IF(grid%label == 'r') THEN
         ALLOCATE ( grid%reg_pt_wt(1:nreg), grid%reg_poly(1:nreg) )
         Call LaGrange_Polynomials(grid,grid%reg_poly)
      ELSE IF( grid%label == 'theta') THEN
!
!        For the angular coordinate, it is necessary to extract the
!        square root behavior of the functions explicitly to get polynomic
!        behavior in the remainder.  This leads to a different approach for
!        even and odd m quantum numbers.
!
!        Compute the Even polynomials.
!
         ALLOCATE ( grid%reg_pt_wt(1:nreg), grid%reg_poly(1:nreg) )
         Call LaGrange_Polynomials(grid,grid%reg_poly)
!
!            Compute the Odd polynomials.
!
         IF (m_max > 0 ) THEN
             ALLOCATE ( grid%reg_poly_odd(1:nreg) )
             Call LaGrange_Polynomials(grid,grid%reg_poly_odd)
         END IF
      END IF
  ELSE IF (keyword == 'spheroidal') THEN
!
!     In the spheroidal coordinate system both the eta and xi coordinates
!     have square root singularities and even and odd m need to be treated
!     differently.
!
!     Compute the Even polynomials.
!  
      ALLOCATE ( grid%reg_pt_wt(1:nreg), grid%reg_poly(1:nreg) )
      Call LaGrange_Polynomials(grid,grid%reg_poly)
!
!      Compute the Odd polynomials.
!
      IF (m_max > 0 ) THEN
          ALLOCATE ( grid%reg_poly_odd(1:nreg) )
          Call LaGrange_Polynomials(grid,grid%reg_poly_odd)
      END IF
  ELSE IF (keyword == 'fourier') THEN
          ALLOCATE ( grid%reg_poly_fourier(1:nreg) )
          Call LaGrange_Polynomials(grid,grid%reg_poly_fourier)
  END IF
1 FORMAT(/,10x,'Calculate the Regional Basis Functions for Coordinate = ',a4)
2 FORMAT(/,10x,'Region = ',i3,2x,'Number of Polynomials = ', i3)
END SUBROUTINE Lobatto_Functions
!***********************************************************************
!***********************************************************************
!deck Coordinate_Factors.f
!***begin prologue     Coordinate_Factors
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute needed prefactors for various
!***                   coordinate systems.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE Coordinate_Factors(grid)
  IMPLICIT NONE
  TYPE(coordinates)                    :: grid
  INTEGER                              :: i
!
  DO i = 1, nreg
    ALLOCATE ( grid%reg_pt_wt(i)%qr_fac( 1:npt(i) ),             &
               grid%reg_pt_wt(i)%inv_qr_fac( 1:npt(i) ),         &
               grid%reg_pt_wt(i)%inv_sqrt_qr_fac( 1:npt(i) ) )
!
!
!     For cartesian coordinates they are all one.
!
    IF (grid%label == 'x' .or. grid%label == 'y' .or. grid%label == 'z' ) THEN  
        grid%reg_pt_wt(i)%qr_fac(:) = one
        grid%reg_pt_wt(i)%inv_qr_fac(:) = one 
        grid%reg_pt_wt(i)%inv_sqrt_qr_fac(:) = one 
    ELSE IF (grid%label == 'r' ) THEN  
!
!     For the radial spherical variable the volume element
!     is r*r.
!
         grid%reg_pt_wt(i)%qr_fac(:) = grid%reg_pt_wt(i)%qr(:) * grid%reg_pt_wt(i)%qr(:) 
         grid%reg_pt_wt(i)%inv_qr_fac(:) = one  / grid%reg_pt_wt(i)%qr_fac(:)
         grid%reg_pt_wt(i)%inv_sqrt_qr_fac(:) = one / grid%reg_pt_wt(i)%qr(:) 
    ELSE IF (grid%label == 'theta'.or.grid%label == 'eta') THEN
!
!     For both the spherical and the spheroidal angular variable
!     we need ( 1 - eta*eta) as the factor where we are using the
!     cosine of the angle as the relevant variable.
!
         grid%reg_pt_wt(i)%qr_fac(:) = one - grid%reg_pt_wt(i)%qr(:) * grid%reg_pt_wt(i)%qr(:)
         grid%reg_pt_wt(i)%inv_qr_fac(:) = one / grid%reg_pt_wt(i)%qr_fac(:)
         grid%reg_pt_wt(i)%inv_sqrt_qr_fac(:) = Sqrt ( grid%reg_pt_wt(i)%inv_qr_fac(:) )
    ELSE IF (grid%label == 'xi' ) THEN
!
!     For the spheroidal radial variable we need ( xi*xi - 1) 
!     as the factor.
!
         grid%reg_pt_wt(i)%qr_fac(:) = grid%reg_pt_wt(i)%qr(:) * grid%reg_pt_wt(i)%qr(:) - one
         grid%reg_pt_wt(i)%inv_qr_fac(:) = one / grid%reg_pt_wt(i)%qr_fac(:)
         grid%reg_pt_wt(i)%inv_sqrt_qr_fac(:) = Sqrt ( grid%reg_pt_wt(i)%inv_qr_fac(:) )
    ELSE IF (grid%label == 'rho') THEN
!
!     For cylindrical coordinates we need rho as the factor.
!
         grid%reg_pt_wt(i)%qr_fac(:) = grid%reg_pt_wt(i)%qr(:)
         grid%reg_pt_wt(i)%inv_qr_fac(:) = one / grid%reg_pt_wt(i)%qr(:)
         grid%reg_pt_wt(i)%inv_sqrt_qr_fac(:) = Sqrt ( grid%reg_pt_wt(i)%inv_qr_fac(:) )
    END IF
  END DO
END SUBROUTINE Coordinate_Factors
!***********************************************************************
!***********************************************************************
!deck FE_DVR_Matrices.f
!***begin prologue     FE_DVR_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute the sector matrix elements
!***                   of the kinetic energy operator.  Note
!***                   that the comments on even and odd m apply
!***                   here as well.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE FE_DVR_Matrices(grid)
  IMPLICIT NONE
  TYPE (coordinates)             :: grid
  INTEGER                        :: i
!
  IF (keyword == 'cartesian') THEN
      ALLOCATE(grid%reg_mat(1:nreg))
      Call Kinetic_Energy(grid, grid%reg_mat, grid%reg_poly)
  ELSE IF(keyword =='spherical') THEN
      IF(grid%label == 'r') THEN
         ALLOCATE ( grid%reg_mat(1:nreg) )
         Call Kinetic_Energy(grid, grid%reg_mat, grid%reg_poly)
      ELSE IF( grid%label == 'theta') THEN
!
!        Compute the Even KE.
!
          pre_factor = - one
         ALLOCATE ( grid%reg_mat(1:nreg) )
         Call Kinetic_Energy(grid, grid%reg_mat, grid%reg_poly)
!
!            Compute the Odd polynomials.
!
         IF (m_max > 0 ) THEN
             ALLOCATE ( grid%reg_mat_odd(1:nreg) )
             Call Kinetic_Energy(grid, grid%reg_mat_odd, grid%reg_poly_odd)
         END IF
      END IF
  ELSE IF(keyword =='cylindrical') THEN
      IF(grid%label == 'rho') THEN
         ALLOCATE ( grid%reg_mat(1:nreg) )
         Call Kinetic_Energy(grid, grid%reg_mat, grid%reg_poly)
      END IF
  ELSE IF (keyword == 'spheroidal') THEN
      ALLOCATE(grid%reg_mat(1:nreg))
      IF (grid%label == 'eta') THEN
          pre_factor = - one
      ELSE IF(grid%label == 'xi') THEN
          pre_factor = one
      END IF
!
!     Compute the even kinetic energy matrix elements.
!
      Call Kinetic_Energy(grid, grid%reg_mat, grid%reg_poly)
      DO i = 1, nreg
         grid%reg_mat(i)%tr(:,:) = R_ab * grid%reg_mat(i)%tr(:,:)
      END DO
      IF ( m_max > 0 ) THEN
           ALLOCATE(grid%reg_mat_odd(1:nreg))
!
!      Compute the odd kinetic energy matrix elements.
!
           Call Kinetic_Energy(grid, grid%reg_mat_odd, grid%reg_poly_odd)
           DO i = 1, nreg
              grid%reg_mat(i)%tr(:,:) = R_ab * grid%reg_mat(i)%tr(:,:)
           END DO
      END IF
  ELSE IF (keyword == 'fourier') THEN
      ALLOCATE( grid%reg_mat_fourier(1))
      Call Kinetic_Energy(grid, grid%reg_mat_fourier)
  ELSE IF (keyword == 'hermite') THEN
      ALLOCATE( grid%reg_mat_hermite(1))
      Call Kinetic_Energy(grid, grid%reg_mat_hermite)
  ELSE IF (keyword == 'laguerre') THEN
      ALLOCATE( grid%reg_mat_laguerre(1))
      Call Kinetic_Energy(grid, grid%reg_mat_laguerre)
  END IF
!
!
END SUBROUTINE FE_DVR_Matrices
!***********************************************************************
!***********************************************************************
           END MODULE FEDVR_Module
!***********************************************************************
!***********************************************************************
