!***********************************************************************
! Two_Electron_FEDVR_Module
!**begin prologue     Two_Electron_FEDVR_Module
!**date written       090219   (yymmdd)
!**revision date               (yymmdd)
!**keywords           DVR, FEDVR
!**
!**author             schneider, b. i.(nsf)
!**source             DVR Library
!**purpose            Calculates two electron integrals in DVR basis
!***                  
!***description       Compute  < rho_ik | 1/r_12 | rho_jl > where
!***                        rho_ij(r) = psi_i(r) * psi_j(r) and 
!***                        psi_i(r_j) = delta_ij/sqrt(w_i)
!***                        The points and weights are those defined
!***                        in the FEDVR quadrature basis. 
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      Two_Electron_FEDVR_Module
!***********************************************************************
!***********************************************************************
                           MODULE Two_Electron_FEDVR_Module
                            USE Data_Module
                            USE FEDVR_Shared
                            USE FEDVR_Derived_Types
                            USE Two_Electron_Shared
                            USE Legendre
!***********************************************************************
!                           Explicit Interfaces
!***********************************************************************
!
!
!
                            INTERFACE Data_Input                  
                       MODULE PROCEDURE Spherical_Data,            &
                                        Spheroidal_Data
                            END INTERFACE Data_Input
!
                           INTERFACE V_ijkl                        
                       MODULE PROCEDURE V_ijkl_Spherical,         &
                                        V_ijkl_Spheroidal  
                            END INTERFACE V_ijkl
!
                           INTERFACE D_LM                       
                       MODULE PROCEDURE D_LM_Spherical_Harmonic,  &
                                        D_LM_FEDVR  
                            END INTERFACE D_LM
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck D_LM_Spherical_Harmonic
!***begin prologue     D_LM_Spherical_Harmonic
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Calculate the D_LM Coefficients
!***description        The D_LM coefficients are products of Wigner 3J coefficients
!***                   (-1)^m * C_3J(l,l',L|0,0,0) * C_3J(l,l',L|-m,m',M)
!***                   and arise from the integation of three spherical harmonics one being
!***                   complex conjugated which leads to the (-1)^m prefactor.
!***references
!***routines called
!***end prologue       D_LM_Spherical_Harmonic
  SUBROUTINE D_LM_Spherical_Harmonic (harmonics)
  IMPLICIT NONE
!  TYPE(angular_matrices)                   :: ang_mat(0:maximum_total_L,0:maximum_total_M)
  TYPE (spherical_harmonics)               :: harmonics
  INTEGER                                  :: l_1
  INTEGER                                  :: l_2
  INTEGER                                  :: m_1
  INTEGER                                  :: m_2  
  INTEGER                                  :: m_fac  
  INTEGER                                  :: l_tot  
  INTEGER                                  :: m_tot
  INTEGER                                  :: l_upper  
  INTEGER                                  :: l_lower
  INTEGER                                  :: count
  INTEGER                                  :: possible
  REAL(idp)                                :: pre_factor
  REAL(idp)                                :: multiplier
  REAL(idp)                                :: F_3J
  REAL(idp)                                :: coef
  REAL(idp)                                :: value
  DO l_1 = 0, l_max
     DO l_2 = 0, l_1
        IF (  ( l_1 + l_2 ) > maximum_total_l ) THEN
           exit
        ELSE
           l_upper = l_1 + l_2
           l_lower = abs(l_1-l_2)
           pre_factor = sqrt ( dfloat ( (l_1 + l_1 + int_one) * ( l_2 + l_2 + int_one ) ) )
           possible = ( l_1 + l_1 + int_one ) * ( l_2 + l_2 + int_one )
           ALLOCATE( ang_mat(l_tot:l_tot,-l_tot:l_tot) )
           DO l_tot = l_lower, l_upper
              coef = pre_factor * F_3J(l_1,int_zero,l_2,int_zero,l_tot,int_zero,.false.) 
              ALLOCATE( ang_mat(l_tot:l_tot,-l_tot:l_tot) )
              DO m_tot = -l_tot, l_tot
                 ALLOCATE( ang_mat(l_tot,m_tot)%D_LM_Coef(1:possible) )
                 m_fac = (-int_one)**(-l_1)
                 multiplier = m_fac * pre_factor
                 count = int_zero
                 DO m_1 = -l_1, l_1
                    DO m_2 = -l_2, l_2
                       IF ( m_tot /= ( m_2-m_1) ) THEN
                            exit
                       ELSE
                          value = multiplier * F_3J(l_1,-m_1,l_2,m_2,l_tot,m_tot,.false.) 
                          IF (value /= zero) THEN
                              count = count + int_one
                              ang_mat(l_tot,m_tot)%D_LM_Coef(count) = value
                              ang_mat(l_tot,m_tot)%D_LM_Index(count,1) = l_1
                              ang_mat(l_tot,m_tot)%D_LM_Index(count,2) = m_1
                              ang_mat(l_tot,m_tot)%D_LM_Index(count,3) = l_2
                              ang_mat(l_tot,m_tot)%D_LM_Index(count,4) = m_2
                          END IF
                       END IF 
                    END DO
                    m_fac = - m_fac
                 END DO
              END DO
           END DO
        END IF
     END DO
  END DO
END SUBROUTINE D_LM_Spherical_Harmonic
!***********************************************************************
!***********************************************************************
!deck D_LM_FEDVR
!***begin prologue     D_LM_FEDVR
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Calculate the D_LM Coefficients
!***description        The D_LM coefficients are products of Wigner 3J coefficients
!***                   (-1)^m * C_3J(l,l',L|0,0,0) * C_3J(l,l',L|-m,m',M)
!***                   and arise from the intergation of thee spherical harmonics one being
!***                   complex conjugated which leads to the (-1)^m prefactor.
!***references
!***routines called
!***end prologue       D_LM_FEDVR
  SUBROUTINE D_LM_FEDVR (fedvr)
  USE P_LM
  IMPLICIT NONE
  TYPE(coordinates), DIMENSION(2)          :: grid
  TYPE(spherical_fedvr)                    :: fedvr
  TYPE(Reg_LM)                             :: R_LM
!  TYPE(angular_matrices)                   :: ang_mat(0:maximum_total_L,0:maximum_total_M)
  INTEGER                                  :: l_1
  INTEGER                                  :: l_2
  INTEGER                                  :: m_1
  INTEGER                                  :: m_2  
  INTEGER                                  :: m_fac  
  INTEGER                                  :: l_tot  
  INTEGER                                  :: m_tot
  INTEGER                                  :: l_upper  
  INTEGER                                  :: l_lower
  INTEGER                                  :: count
  INTEGER                                  :: possible
  INTEGER                                  :: len
  REAL(idp)                                :: pre_factor
  REAL(idp)                                :: multiplier
  REAL(idp)                                :: F_3J
  REAL(idp)                                :: coef
  REAL(idp)                                :: value
  Call pakstr(FEDVR_File,len)
  write(iout,*) 'Opening Data File = '//FEDVR_File//' as old'
  Call IOsys('open '//FEDVR_File//' as old',0,0,0,FEDVR_File)
  Call IOsys('read real "number of physical points" from '//FEDVR_File,1,        &
              number_of_angular_points,0,' ')
  write(iout,1) number_of_angular_points
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
                                          number_of_angular_points,1,            &
                                          iout,'e')
  ALLOCATE( Factorial(0:maximum_total_l + maximum_total_m) )
  Call N_Factorial 
  ALLOCATE(Leg%R_LM%F_x(0:maximum_total_l,0:maximum_total_m,1:number_of_angular_points))
  n_points=number_of_angular_points
  Call Legendre(R_LM, grid(2)%grid_points, n_points)
1 Format(/,15x,'number of angular points from disk       = ',i5)
END SUBROUTINE D_LM_FEDVR
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
  IF (representation == 'fedvr') THEN
      FEDVR_File = 'spherical'//'_'//'theta'
      Call D_LM(fedvr)
  ELSE IF (representation == 'spherical harmonics') THEN
      Call D_LM(harmonics)
  END IF
  ALLOCATE(dvr_mat(0:l_max))
  Call V_ijkl_Spherical (grid(1), name_spherical)
  write(iout,*) 'Opening Data File = Spherical_2_Electron_Integrals as new'
  Call IOsys('open FEDVR_Two_Electron_Integral_File as new',0,0,0,                  &
             'Spherical_2_Electron_Integrals')
  Call IOsys('rewind all on '//FEDVR_File//' read-and-write',0,0,0,' ')
1 Format(/,15x,'maximum orbital l from disk       = ',i4,/,15x,                      &
               'maximum orbital m from disk       = ',i4,/,15x,                      &
               'number of radial points from disk = ',i5,/,15x,                      &
               'box size from disk                = ',f15.8)

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
!deck V_ijkl_Spherical
!***begin prologue     V_ijkl_Spherical
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            V_ijkl_Spherical
!***description        The radial parts of the two-electron integrals
!***                   in spherical coordinates are computed for all
!***                   required L values.  They do not explicitly depend on M.
!***references
!***routines called
!***end prologue       V_ijkl_Spherical
  SUBROUTINE V_ijkl_Spherical (grid, name_spherical)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE(spherical)                          :: name_spherical
  REAL(idp)                                :: R_factor  
  REAL(idp), DIMENSION(:), ALLOCATABLE     :: factr
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: l_factor  
  INTEGER                                  :: l_two
  INTEGER                                  :: info  
  INTEGER                                  :: l_1  
  INTEGER                                  :: l_2  
  INTEGER                                  :: l_3  
  INTEGER                                  :: lm  
  INTEGER                                  :: count
  CHARACTER(LEN=3)                         :: itoc
  CHARACTER(LEN=80)                        :: title
  REAL(idp)                                :: value
!
!       V_ijkl = delta_ik delta_jl ( r_i )^(L+2)*( r_j )^(L+2) /( r_n)^(2L+1) 
!                                          - 
!                         (2L+1) * ( r_i )^2 * ( T_ij )^-1 /sqrt(w_i*w_j}
!
!
!      Compute the inverse from the LU decomposition.  The last point
!      is excluded since the solution of the Poisson equation we are
!      computing is zero at the last point.  We will add a solution of the
!      homogeneous equation to take care of the boundary conditions later. 
!      
!         
  ALLOCATE( dvr_mat(0)%ipvt( n_final ), dvr_mat(0)%lower( n_tri ),                  &
            dvr_mat(0)%Inverse(n_final,n_final),                                    &
            dvr_mat(0)%tr(physical_points,physical_points), factr( n_total ) )
  factr(1:n_total) = grid%grid_points(1:n_total) / sqrt ( grid%grid_weights(1:n_total) ) 
  keyword = 'spherical'
  l_1=lenth(keyword)
  grid%label = 'r'
  l_2=lenth(grid%label)
!
!    Set up the matrix to be used to solve the linear equations
!
  DO lm = 0, maximum_total_L
     title=keyword(1:l_1)//'_'//grid%label(1:l_2)//'_'//itoc(lm)
     l_3=lenth(title)
     write(iout,*) 'Reading Global Nabla = ',title(1:l_3)//' from disk'
     Call IOsys('read real "'//title(1:l_3)//' nabla" from '//FEDVR_File,             &
                 physical_points*physical_points,dvr_mat(0)%tr,0,' ')
     count = 0
     DO i = 1, n_final
        DO j = 1, i
           count = count + 1
           dvr_mat(0)%lower(count) = dvr_mat(lm)%tr(i,j)
        END DO
     END DO
!
!    Factor the matrix
!
     Call DSPTRF ('u',n_final,dvr_mat(0)%lower,dvr_mat(0)%ipvt,info)

!
!    Compute the inverse
!
     dvr_mat(0)%Inverse(:,:) = zero
     DO i = 1, n_final
        dvr_mat(0)%Inverse(i,i) = one
     END DO
!
     Call DSPTRS('u',n_final,n_final,dvr_mat(0)%lower,dvr_mat(0)%ipvt,             &
                                     dvr_mat(0)%Inverse,n_total,info)
!
!         
     DO i = 1, n_final
        dvr_mat(0)%Inverse(1:n_final,i) = factr(1:n_final) * dvr_mat(0)%Inverse(1:n_final,i) * factr(i) 
     END DO
     l_factor = - ( lm + lm + int_one )
     dvr_mat(0)%Inverse(:,:) = l_factor * dvr_mat(0)%Inverse(:,:)
     l_two = lm + int_two     
     factr(1:n_total) = grid%grid_points(1:n_total) ** l_two
!
!    Add the solution to the homogeneous equation
!                                                                           
     ALLOCATE( dvr_mat(lm)%Q(n_total,n_total) )
     write(iout,*) 'Begin Calculation of Basic Radial Two-Electron Integrals L = '//itoc(lm)
     dvr_mat(lm)%Q(:,:) = zero
     dvr_mat(lm)%Q(1:n_final,1:n_final) = dvr_mat(0)%Inverse(1:n_final,1:n_final)
     R_factor = one / R_max**l_factor
     DO i = 1, n_total
        dvr_mat(lm)%Q(1:n_total,i) = dvr_mat(lm)%Q(1:n_total,i) + R_factor * factr(1:n_total) * factr(i)
     END DO
     write(iout,*) 'End Calculation of Basic Radial Two-Electron Integrals L = '//itoc(lm)
     title='Q Function for Two Electron Integrals_'//itoc(lm)
     Call Prntfm(title,dvr_mat(lm)%Q,n_total, n_total, n_total, n_total, iout,'e')
  END DO
  DEALLOCATE( dvr_mat(0)%ipvt, dvr_mat(0)%lower, dvr_mat(0)%Inverse, dvr_mat(0)%tr )
END SUBROUTINE V_ijkl_Spherical
!***********************************************************************
!***********************************************************************
!deck V_ijkl_Spheroidal
!***begin prologue     V_ijkl_Spheroidal
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Q function for computation of two-electron integrals.
!***
!***references
!***routines called
!***end prologue       V_ijkl_Spheroidal
  SUBROUTINE V_ijkl_Spheroidal(grid, val, name_spheroidal)
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
END SUBROUTINE V_ijkl_Spheroidal

!***********************************************************************
!***********************************************************************
           END MODULE Two_Electron_FEDVR_Module
!***********************************************************************
!***********************************************************************
