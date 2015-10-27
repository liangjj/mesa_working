!***********************************************************************
! Two_Electron_Integral_Setup_Module
!**begin prologue     Two_Electron_Integral_Setup_Module
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             
!**purpose            Read in FEDVR Data
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      Two_Electron_Integral_Setup_Module
!***********************************************************************
!***********************************************************************
                           MODULE Two_Electron_Integral_Setup_Module
                           USE Data_Module
                           USE FEDVR_Shared
                           USE FEDVR_Derived_Types
                           USE Two_Electron_Shared
!                           USE Two_Electron_FEDVR_Integrals
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                            INTERFACE Data_Input                       
                       MODULE PROCEDURE Spherical_Data,                        &
                                        Spheroidal_Data
                            END INTERFACE Data_Input
!
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Setup_2e
!***begin prologue     Setup_2e
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            setup for two electron integral calculation
!***description        
!***references
!***routines called    
!***end prologue       Setup_2e
  SUBROUTINE Setup_2e  
  IMPLICIT NONE
  REAL(idp)                   :: fpkey
  LOGICAL                     :: dollar
  LOGICAL                     :: logkey
  INTEGER                     :: intkey  
  CHARACTER(LEN=80)           :: chrkey
  CHARACTER(LEN=8)            :: itoc
!
  IF ( dollar('$Two_Electron_Integrals',card,cpass,inp) )THEN
       keyword = chrkey(card,'coordinate_system','spherical',' ')
       len=lenth(keyword)
       IF (keyword(1:len) /= 'spherical' .or. keyword(1:len) /= 'spheroidal') THEN
           Call lnkerr('Quit.  Coordinate System not Spherical or Spheroidal')           
       END IF
       write(iout,1) keyword(1:len)
       representation = chrkey(card,'representation','spherical_harmonics',' ')
       maximum_orbital_l = intkey(card,'maximum_orbital_l',0,' ')
       maximum_orbital_m = intkey(card,'maximum_orbital_m',maximum_orbital_l,' ')
       minimum_orbital_m = - maximum_orbital_l 
       maximum_total_L = intkey(card,'maximum_total_L',maximum_orbital_l + maximum_orbital_l,' ') 
       maximum_total_M = intkey(card,'maximum_total_M',maximum_total_L,' ')  
       minimum_total_M = - maximum_total_M
       ALLOCATE(reg_grid(2))
  END IF
  write(iout,1) maximum_orbital_l, maximum_orbital_m, maximum_total_L, maximum_total_M
  IF (keyword(1:len) == 'spherical') THEN
      Call Data_Input(reg_grid,reg_grid(1)%name_spherical)
  ELSE IF (keyword(1:len) == 'spheroidal') THEN
      Call Data_Input(reg_grid,reg_grid(1)%name_spheroidal)
  END IF
1 FORMAT(/,25x,'Coordinate System = ',a32)
2 Format(/,15x,'maximum orbital l read in = ',i4,/,15x,                               &
               'maximum orbital m read in = ',i4,/,15x,                               &
               'maximum total L read in   = ',i4,/,15x,                               &
               'maximum total M read in   = ',i4 )
  END SUBROUTINE Setup_2e
!***********************************************************************
!***********************************************************************
!deck Spherical_Data
!***begin prologue     Spherical_Data
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine 
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       Spherical_Data
  SUBROUTINE Spherical_Data(grid,name_spherical)
  IMPLICIT NONE
  TYPE(coordinates), DIMENSION(2)  :: grid
  TYPE(spherical)                  :: name_spherical
  TYPE(spherical_harmonics)        :: harmonics
  TYPE(spherical_fedvr)            :: fedvr
  CHARACTER(LEN=3)                 :: itoc
  INTEGER                          :: lm
  INTEGER                          :: len_1
  INTEGER                          :: l_1
  INTEGER                          :: l_2
!
  FEDVR_File = 'spherical'//'_'//'r'
  Call pakstr(FEDVR_File,len_1)
  write(iout,*) 'Opening Data File = '//FEDVR_File//' as old'
  Call IOsys('open '//FEDVR_File//' as old',0,0,0,FEDVR_File)
  Call IOsys('read integer "L Max" from '//FEDVR_File,1,l_max,0,' ')
  Call IOsys('read integer "M Max" from '//FEDVR_File,1,m_max,0,' ')
  Call IOsys('read real "number of physical points" from '//FEDVR_File,1,            &
              number_of_radial_points,0,' ')
  Call IOsys('read real "drop right" from '//FEDVR_File,1,                           &
              drop(2),0,' ')
  Call IOsys('read real "last grid point" from '//FEDVR_File,1,R_max,0,' ')
  n_final = number_of_radial_points
  n_total = number_of_radial_points
  IF(.not.drop(2)) THEN
      n_final = number_of_radial_points - int_one
  ELSE
      n_total = number_of_radial_points + int_one
  END IF
  n_tri = n_final * (n_final + int_one ) / int_two
  ALLOCATE(grid(1)%grid_points(number_of_radial_points),                             &
           grid(1)%grid_weights(number_of_radial_points) ) 
  Call IOsys('read real "physical grid points" from '//FEDVR_File,                   &
              number_of_radial_points,grid(1)%grid_points,0,' ')
  Call IOsys('read real "physical grid weights" from '//FEDVR_File,                  &
              number_of_radial_points,grid(1)%grid_weights,0,' ')
  Call IOsys('close '//FEDVR_File,0,0,0,' ')
  write(iout,1) l_max, m_max, number_of_radial_points, R_max
  title = 'Radial Points'
  Call prntfmn(title,grid(1)%grid_points,number_of_radial_points,1,                  &
                                         number_of_radial_points,1,                  &
                                         iout,'e')
  title = 'Radial Weights'
  Call prntfmn(title,grid(1)%grid_weights,number_of_radial_points,1,                 &
                                          number_of_radial_points,1,                 &
                                          iout,'e')
  IF (representation /= 'spherical_harmonics') THEN
      FEDVR_File = 'spherical'//'_'//'theta'
      Call pakstr(FEDVR_File,len_1)
      write(iout,*) 'Opening Data File = '//FEDVR_File//' as old'
      Call IOsys('open '//FEDVR_File//' as old',0,0,0,FEDVR_File)
      Call IOsys('read real "number of physical points" from '//FEDVR_File,1,        &
                  number_of_angular_points,0,' ')
      write(iout,2) number_of_angular_points
      ALLOCATE(grid(2)%grid_points(number_of_angular_points),                        &
               grid(2)%grid_weights(number_of_angular_points) ) 
      Call IOsys('read real "physical grid points" from '//FEDVR_File,               &
                  number_of_angular_points,grid(2)%grid_points,0,' ')
      Call IOsys('read real "physical grid weights" from '//FEDVR_File,              &
                  number_of_angular_points,grid(2)%grid_weights,0,' ')
      Call IOsys('close '//FEDVR_File,0,0,0,' ')
      title = 'Angular Points'
      Call prntfmn(title,grid(2)%grid_points,number_of_angular_points,1,             &
                                             number_of_angular_points,1,             &
                                             iout,'e')
      title = 'Angular Weights'
      Call prntfmn(title,grid(2)%grid_weights,number_of_angular_points,1,            &
                                              number_of_radial_points,1,             &
                                              iout,'e')
  END IF
  ALLOCATE(dvr_mat(0:l_max))
  Call V_ijkl_Spherical (grid(1), name_spherical)
  IF (representation == 'spherical_harmonics') THEN
      Call D_LM(harmonics)
  ELSE
      Call D_LM(fedvr)
  END IF
  write(iout,*) 'Opening Data File = Spherical_2_Electron_Integrals as new'
  Call IOsys('open FEDVR_Two_Electron_Integral_File as new',0,0,0,                  &
             'Spherical_2_Electron_Integrals')
  Call IOsys('rewind all on '//FEDVR_File//' read-and-write',0,0,0,' ')
1 Format(/,15x,'maximum orbital l from disk       = ',i4,/,15x,                      &
               'maximum orbital m from disk       = ',i4,/,15x,                      &
               'number of radial points from disk = ',i5,/,15x,                      &
               'box size from disk                = ',f15.8)
2 Format(/,15x,'number of angular points from disk       = ',i5)
  END SUBROUTINE Spherical_Data
!***********************************************************************
!***********************************************************************
!deck Spheroidal_Data
!***begin prologue     Spheroidal_Data
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine 
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       Spheroidal_Data
  SUBROUTINE Spheroidal_Data(grid,name_spheroidal)
  IMPLICIT NONE
  TYPE(coordinates), DIMENSION(2)  :: grid
  TYPE(spheroidal)                 :: name_spheroidal
  INTEGER                          :: len_1
!
  FEDVR_File = 'spheroidal'//'_'//'eta'
  Call pakstr(FEDVR_File,len_1)
  write(iout,*) 'Opening Data File = '//FEDVR_File//' as old'
  Call IOsys('open '//FEDVR_File//' as old',0,0,0,FEDVR_File)
  Call IOsys('read integer "L Max" from '//FEDVR_File,1,l_max,0,' ')
  Call IOsys('read integer "M Max" from '//FEDVR_File,1,m_max,0,' ')
  Call IOsys('read real "number of physical points" from '//FEDVR_File,1,            &
              number_of_radial_points,0,' ')
  Call IOsys('read real "drop right" from '//FEDVR_File,1,                           &
              drop(2),0,' ')
  Call IOsys('read real "last grid point" from '//FEDVR_File,1,R_max,0,' ')
  n_final = number_of_radial_points
  n_total = number_of_radial_points
  IF(.not.drop(2)) THEN
      n_final = number_of_radial_points - int_one
  ELSE
      n_total = number_of_radial_points + int_one
  END IF
  n_tri = n_final * (n_final + int_one ) / int_two
  ALLOCATE(grid(1)%grid_points(number_of_radial_points),                             &
           grid(1)%grid_weights(number_of_radial_points) ) 
  Call IOsys('read real "physical grid points" from '//FEDVR_File,                   &
              number_of_radial_points,grid(1)%grid_points,0,' ')
  Call IOsys('read real "physical grid weights" from '//FEDVR_File,                  &
              number_of_radial_points,grid(1)%grid_weights,0,' ')
  Call IOsys('close '//FEDVR_File,0,0,0,' ')
  FEDVR_File = 'spherical'//'_'//'xi'
  Call pakstr(FEDVR_File,len_1)
  write(iout,*) 'Opening Data File = '//FEDVR_File//' as old'
  Call IOsys('open '//FEDVR_File//' as old',0,0,0,FEDVR_File)
  Call IOsys('read real "number of physical points" from '//FEDVR_File,1,            &
              number_of_angular_points,0,' ')
  ALLOCATE(grid(2)%grid_points(number_of_angular_points),                            &
           grid(2)%grid_weights(number_of_angular_points) ) 
  Call IOsys('read real "physical grid points" from '//FEDVR_File,                   &
              number_of_angular_points,grid(2)%grid_points,0,' ')
  Call IOsys('read real "physical grid weights" from '//FEDVR_File,                  &
              number_of_radial_points,grid(2)%grid_weights,0,' ')
  Call IOsys('close '//FEDVR_File,0,0,0,' ')
  END SUBROUTINE Spheroidal_Data
!***********************************************************************
!***********************************************************************
           END MODULE Two_Electron_Integral_Setup_Module
!***********************************************************************
!***********************************************************************
