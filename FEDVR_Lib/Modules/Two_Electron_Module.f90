!***********************************************************************
! Two_Electron_Module
!**begin prologue     Two_Electron_Module
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
!***end prologue      Two_Electron_Module
!***********************************************************************
!***********************************************************************
                           MODULE Two_Electron_Module
                           USE Data_Module
                           USE FEDVR_Shared
                           USE FEDVR_Derived_Types
                     INTEGER                      :: len_1
!***********************************************************************
!                           Explicit Interfaces
!***********************************************************************
!
!
                            INTERFACE Two_Electron_Q                       
                       MODULE PROCEDURE Two_Electron_Spherical_Q,         &
                                        Two_Electron_Spheroidal_Q  
                            END INTERFACE Two_Electron_Q  
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Two_Electron_Driver
!***begin prologue     Two_Electron_Driver
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Solve Poisson equation
!***                   
!***                   
!***                   
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Two_Electron_Driver
!
  SUBROUTINE Two_Electron_Driver(grid,label)
  IMPLICIT NONE
  TYPE (coordinates)                        :: grid
  CHARACTER(LEN=*)                          :: label
  INTEGER                                   :: i
  INTEGER                                   :: lm
  LOGICAL                                   :: dollar
!  
!
  IF ( dollar('$two_electron_'//grid%label(1:len_1),card,cpass,inp) ) THEN
       two_electron=.true. 
  ELSE
       two_electron=.false. 
  END IF
  IF ( two_electron ) THEN 
       n_total=physical_points
      IF(.not.drop(2)) THEN
          n_total=n_total - 1
      END IF
      ALLOCATE( dvr_mat(0)%ipvt( n_total ), dvr_mat(0)%lower( n_total*(n_total+1)/2 ),   &
                dvr_mat(0)%work( 5*n_final ) )
  END IF
  DO lm = 0, size
     ALLOCATE( dvr_mat(lm)%tr( physical_points, physical_points),                        &          
               dvr_mat(lm)%Q( physical_points, physical_points) )                                   
     Call Solver(grid,lm)
     IF(keyword == 'spherical') THEN
         Call Two_Electron_Q(grid,lm,grid%name_spherical)
     ELSE IF(keyword == 'spheroidal') THEN
         Call Two_Electron_Q(grid,lm,grid%name_spheroidal)
     END IF
     DEALLOCATE(dvr_mat(lm)%ham)
     len=lenth(keyword)
     len_1=lenth(label)
     title=keyword(1:len)//'_'//label(1:len_1)
     len=lenth(title)
     Call IOsys('write real "'//title(1:len)//' Q-matrix" to '//FEDVR_File,                &
                 physical_points*number_of_right_hand_sides,dvr_mat(0)%Q,0,' ')
  END DO
END SUBROUTINE Two_Electron_Driver
!***********************************************************************
!***********************************************************************
!deck Two_Electron_Spherical_Q
!***begin prologue     Two_Electron_Spherical_Q
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Q function for computation of two-electron integrals
!***
!***references
!***routines called
!***end prologue       Two_Electron_Spherical_Q
  SUBROUTINE Two_Electron_Spherical_Q (grid, val, name_spherical)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE(spherical)                          :: name_spherical
  REAL(idp)                                :: R_factor  
  REAL(idp)                                :: fac_1
  REAL(idp)                                :: fac_2  
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: val
  INTEGER                                  :: l_factor  
  CHARACTER(LEN=3)                         :: itoc
  CHARACTER(LEN=80)                        :: title
!
! Since we are interested in the solution where the right hand side is the unit matrix, we need to make
! sure the r**2 in the volume is handled properly.  This ensure that to be the case.
  DO i = 1, physical_points
     DO j = 1, physical_points
        dvr_mat(val)%tr(i,j) = dvr_mat(val)%tr(i,j) / ( grid%grid_points(i) *  grid%grid_points(j) )  
     END DO
  END DO
!
!       Q^L(i,j) = ( r_i )^(L+2)*( r_j )^L /( r_n)^(2L+1) - (2L+1) * ( r_i )^2 * ( T_ji )^-1 /sqrt(w_i*w_j}
!
  l_factor = ( val + val + int_one )
  R_factor = one / R_max**l_factor
  DO i = 1, physical_points
     fac_1 = ( grid%grid_points(i) ) ** ( val + int_two ) * R_factor 
     fac_2 = l_factor * ( grid%grid_points(i) ) ** int_two
     DO j = 1, physical_points
        dvr_mat(val)%Q(i,j) = ( grid%grid_points(j) ** val ) * fac_1 - fac_2 * dvr_mat(val)%tr(j,i)    &
                                                       /                                               &
                                ( sqrt ( grid%grid_weights(i) * grid%grid_weights(j) ) )
     END DO
  END DO
  title='Q Function for Two Electron Integrals_'//itoc(val)
  Call Prntfm(title,dvr_mat(val)%Q,physical_points, physical_points, physical_points,                  &
              physical_points, iout,'e')
END SUBROUTINE Two_Electron_Spherical_Q
!***********************************************************************
!***********************************************************************
!deck Two_Electron_Spheroidal_Q
!***begin prologue     Two_Electron_Spheroidal_Q
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Q function for computation of two-electron integrals.
!***
!***references
!***routines called
!***end prologue       Two_Electron_Spheroidal_Q
  SUBROUTINE Two_Electron_Spheroidal_Q(grid, val, name_spheroidal)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE(spheroidal)                         :: name_spheroidal
  REAL(idp)                                :: scale  
  INTEGER                                  :: info
  INTEGER                                  :: count
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: val
  CHARACTER(LEN=3)                         :: itoc
  CHARACTER(LEN=80)                        :: title
!
END SUBROUTINE Two_Electron_Spheroidal_Q

!***********************************************************************
!***********************************************************************
           END MODULE Two_Electron_Module
!***********************************************************************
!***********************************************************************
