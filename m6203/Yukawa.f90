!deck Yukawa
!***begin prologue     Yukawa
!***date written       140601   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Yukawa
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            Read in radial and angular information needed to compute
!***                   the radial points and weights.
!***description  
!***             
!***             
!***references
!***routines called    
!***end prologue       Yukawa
!********************************************************************************
!********************************************************************************
                        MODULE Yukawa
  USE Data
  USE Grid_Defined_Types
  USE Matrix_Print
  IMPLICIT NONE
!
!********************************************************************************
!********************************************************************************
                       Contains
!***********************************************************************
!***********************************************************************
!deck Yukawa_Potential
!***begin prologue     Yukawa_Potential
!***date written       930502   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           yukawa, link m6200
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            calculate multi-center potential'
!***description        Each atomic cneter has a yukawa potential and it is tabulated
!***                   about that atomic center. 
!***references         
!
!***routines called
!***end prologue       Yukawa_Potential
  subroutine Yukawa_Potential (atom,shl,center)
  IMPLICIT NONE
  TYPE(ATOMS)                                        :: atom
  TYPE(ATOMS), DIMENSION(:)                          :: center
  TYPE(SHELLS)                                       :: shl
  REAL(idp)                                          :: sdot
  REAL(idp)                                          :: dist
  INTEGER                                            :: i
  INTEGER                                            :: nc
  TYPE(REAL_MATRIX)                                  :: type_real_matrix
  TYPE(REAL_VECTOR)                                  :: type_real_vector
!
! For each point in the shell compute the total potential due to the sum of Yukawa potentials from
! all the other atoms.
!
  ALLOCATE(shl%yukawa(1:shl%n_3d) )  
  DO i = 1, shl%n_3d
     DO nc = 1, ncent
        dist = sqrt (                                            &
                    ( shl%xyz_grid(i,1)-center(nc)%cen(1) )  *   &
                    ( shl%xyz_grid(i,1)-center(nc)%cen(1) )  +   &
                    ( shl%xyz_grid(i,2)-center(nc)%cen(2) )  *   &
                    ( shl%xyz_grid(i,2)-center(nc)%cen(2) )  +   &
                    ( shl%xyz_grid(i,3)-center(nc)%cen(3) )  *   &
                    ( shl%xyz_grid(i,3)-center(nc)%cen(3) )  )
        shl%yukawa(i) = shl%yukawa(i) + exp( - center(nc)%eta * dist )  / dist        
     END DO
  END DO
  atom%yukawa_integral = atom%yukawa_integral                &
                                 +                           &
                         sdot(shl%n_3d,shl%yukawa,1,shl%integration_weight,1)
  END Subroutine Yukawa_Potential
!***********************************************************************
!***********************************************************************
END MODULE Yukawa
!***********************************************************************
!***********************************************************************
