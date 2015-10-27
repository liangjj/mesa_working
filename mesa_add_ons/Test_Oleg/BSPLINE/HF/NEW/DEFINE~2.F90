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

 
