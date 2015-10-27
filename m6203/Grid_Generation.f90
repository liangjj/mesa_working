!deck Grid_Generation
!***begin prologue     Grid_Generation
!***date written       140601   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Grid_Generation
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            Read in radial and angular information needed to compute
!***                   the radial points and weights.
!***description  
!***             
!***             
!***references
!***routines called    
!***end prologue       Grid_Generation
!********************************************************************************
!********************************************************************************
                        MODULE Grid_Generation
  USE Data
  USE Grid_Defined_Types
  USE Matrix_Print
  USE Lebedev_Quadrature
  USE Gauss_Quadrature
  USE Renormalization
  IMPLICIT NONE
!
!********************************************************************************
!********************************************************************************
!
                            INTERFACE Quadrature_Grids
                       MODULE PROCEDURE Lebedev_Grid,               &
                                        Angular_Grid,               &  
                                        Radial_Grid  
                            END INTERFACE Quadrature_Grids
!
                            INTERFACE Total_Points
                       MODULE PROCEDURE Lebedev_Points,             &
                                        Gauss_Points
                            END INTERFACE Total_Points
!
                            INTERFACE Set_Grid_Variables
                       MODULE PROCEDURE Set_Angular_Variables,      &  
                                        Set_Radial_Variables  
                            END INTERFACE Set_Grid_Variables
!
                            INTERFACE Fix_BC
                       MODULE PROCEDURE Fix_Angular_BC,             &  
                                        Fix_Radial_BC  
                            END INTERFACE Fix_BC
!
                            INTERFACE Atomic_Grid
                       MODULE PROCEDURE Three_d_Lebedev_Grid,           &  
                                        Three_d_Gauss_Grid  
                            END INTERFACE Atomic_Grid
!
!********************************************************************************
                                 Contains
!********************************************************************************
!********************************************************************************
!deck satshl.f
!***begin prologue     satshl
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           satshl, link m6200
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            read in and write out atomic quadrature information
!***                   This includes points, weights and FEDVR information on functions and derivatives.
!***references

!***routines called
!***end prologue       satshl
  subroutine satshl(atom,center)
  IMPLICIT NONE
  TYPE(ATOMS)                            :: atom
  TYPE(ATOMS), DIMENSION(:)              :: center
  INTEGER                                :: lenth 
  INTEGER                                :: intkey
  INTEGER                                :: precision
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: ii
  INTEGER                                :: i_sign
  INTEGER                                :: j_sign
  INTEGER                                :: rule
  INTEGER                                :: diff
  INTEGER                                :: len
  INTEGER                                :: ang_len
  INTEGER                                :: rad_len
  INTEGER                                :: ang_loop
  CHARACTER (LEN=3)                      :: itoc
  CHARACTER (LEN=3)                      :: yes_no
  CHARACTER (LEN=80)                     :: chrkey
  LOGICAL                                :: logkey
  LOGICAL                                :: dollar
  REAL(idp)                              :: start
  REAL(idp)                              :: end
  REAL(idp)                              :: increment
  REAL(idp)                              :: fpkey
!
!
!                   
!
  ALLOCATE(atom%local_time(1:3),atom%total_time(1:1))
  call pakstr(ang_grid,ang_len)     
  IF( ang_grid(1:ang_len) == 'yes' ) THEN
      nonsep=logkey(card,'angular_quadrature=lebedev',.false.,' ')
      write(iout,*)
      write(iout,*)
      write(iout,*) '                             Begin Angular Coordinates'
      IF (nonsep == .true.) THEN
          t_start=secnds(0.0)
          write(iout,1)
          atom%coord_label='lebedev'
          yes_no = 'yes'
          IF (iosys_on == .true.) THEN
              CALL iosys('write character "non separable quadrature" to grid',0,0,0,'yes')
          END IF
          call pakstr(atom%coord_label,len)     
          ang_loop=atom%n_shl
          IF (fixed_angular_quadrature == .true. ) THEN
              ang_loop=1
          END IF
          DO i = 1, ang_loop
             str="$"//atom%coord_label(1:len)//'_shell_'//itoc(i)
             IF ( dollar(str,card,cpass,inp) ) THEN  
                  write(iout,2) i
                  Call Quadrature_Grids(atom%shl(i)%ang%leb_ang)
                  atom%shl(i)%ang%leb_ang%n_pts=atom%shl(i)%ang%leb_ang%nang
             ELSE
                  Call lnkerr('error in processing card')
             END IF
          END DO
          atom%local_time(1) = secnds(0.0) - t_start
          write(iout,*) 'Time for Lebedev angular grid generation = ',atom%local_time(1)
      ELSE
          t_start=secnds(0.0)
          write(iout,4)
          yes_no = 'no'
          IF (iosys_on == .true.) THEN
              CALL iosys('write character "non separable quadrature" to grid',0,0,0,'no')
          END IF
          ang_loop=atom%n_shl
          IF (fixed_angular_quadrature == .true. ) THEN
              ang_loop=1
          END IF
          atom%coord_label='theta'     
          call pakstr(atom%coord_label,len)     
          str="$"//atom%coord_label(1:len)//'_boundary_conditions'
          IF ( dollar(str(1:len),card,cpass,inp) ) THEN  
               Call Set_BC(atom%coord_label,atom%theta_ang%bc)
          ELSE
               Call lnkerr('error in processing card')
          END IF
          atom%coord_label='phi'     
          call pakstr(atom%coord_label,len)     
          str="$"//atom%coord_label(1:len)//'_boundary_conditions'
          IF ( dollar(str(1:len),card,cpass,inp) ) THEN  
               Call Set_BC(atom%coord_label,atom%phi_ang%bc)
          ELSE
               Call lnkerr('error in processing card')
          END IF
          DO i = 1, ang_loop      ! Loop over the shells
             atom%coord_label='theta'     
             call pakstr(atom%coord_label,len)     
             write(iout,5) i
             str="$"//atom%coord_label(1:len)//'_shell_'//itoc(i)
             call pakstr(str,len)     
             IF ( dollar(str(1:len),card,cpass,inp) ) THEN  
                  Call Set_Grid_Variables(atom,atom%shl(i)%ang,atom%shl(i)%ang%theta_ang%reg)
                  atom%shl(i)%ang%theta_ang%n_reg=atom%shl(i)%ang%n_reg
                  Call Fix_BC(atom%shl(i)%ang,atom%shl(i)%ang%theta_ang%reg,atom%theta_ang%bc)
                  atom%shl(i)%ang%theta_ang%n_phy=atom%shl(i)%ang%n_phy
                  Call Quadrature_Grids(atom,atom%shl(i)%ang,atom%shl(i)%ang%theta_ang%reg)
                  atom%shl(i)%ang%theta_ang%n_pts = atom%shl(i)%ang%n_pts
             ELSE
                  Call lnkerr('error in processing card')
             END IF
             atom%coord_label='phi'     
             call pakstr(atom%coord_label,len)     
             write(iout,6) i
             str="$"//atom%coord_label(1:len)//'_shell_'//itoc(i)
             call pakstr(str,len)     
             IF ( dollar(str(1:len),card,cpass,inp) ) THEN  
                  Call Set_Grid_Variables(atom,atom%shl(i)%ang,atom%shl(i)%ang%phi_ang%reg)
                  atom%shl(i)%ang%phi_ang%n_reg=atom%shl(i)%ang%n_reg
                  Call Fix_BC(atom%shl(i)%ang,atom%shl(i)%ang%phi_ang%reg,atom%phi_ang%bc)
                  atom%shl(i)%ang%phi_ang%n_phy=atom%shl(i)%ang%n_phy
                  Call Quadrature_Grids(atom,atom%shl(i)%ang,atom%shl(i)%ang%phi_ang%reg)
                  atom%shl(i)%ang%n_pts = atom%shl(i)%ang%n_pts * atom%shl(i)%ang%theta_ang%n_pts
             ELSE
                  Call lnkerr('bad coordinate card')
             END IF
             write(iout,*)
             write(iout,7) i
             write(iout,*)
          END DO
          atom%local_time(1) = secnds(0.0) - t_start
          write(iout,*) 'Time for Gauss angular grid generation = ',atom%local_time(1)
      END IF
  END IF 
  atom%total_time(1) = atom%local_time(1)
  call pakstr(rad_grid,rad_len)     
  IF( rad_grid(1:rad_len) == 'yes' ) THEN
      t_start=secnds(0.0)
      atom%coord_label='radial'     
      call pakstr(atom%coord_label,len)     
      str="$"//atom%coord_label(1:len)//'_boundary_conditions'
      IF ( dollar(str(1:len),card,cpass,inp) ) THEN  
           Call Set_BC(atom%coord_label,atom%rad%bc)
      ELSE
           Call lnkerr('error in processing card')
      END IF
      Call Set_Grid_Variables(atom,atom%rad,atom%rad%reg)
!
!  Fix up the first and the last shell.
!
      Call Fix_BC(atom,atom%rad%reg,atom%rad%bc)
      Call Quadrature_Grids( atom, atom%rad, atom%rad%reg)
!
      write(iout,*) '              End Radial Coordinate'
      atom%local_time(2) = secnds(0.0) - t_start
      write(iout,*) 'Time for Gauss radial grid generation = ',atom%local_time(2)
      atom%total_time(1) = atom%total_time(1) + atom%local_time(2)
  END IF
  IF ( make_3d_grid /= .true.) THEN
       return
  END IF      
!
!  Now for each shell create the 3D grid
!
  IF( ang_grid(1:ang_len) == 'yes' .and. rad_grid(1:rad_len) == 'yes') THEN
      t_start=secnds(0.0)
      IF (iosys_on == .true.) THEN
          CALL iosys('read character "non separable quadrature" from grid',-1,0,0,yes_no)
      END IF
      atom%vol=0.d0
      atom%yukawa_integral = 0.d0
      yukawa_integral = 0.d0
      atom%n_pts = 0
      IF (yes_no == 'yes') THEN
!         Lebedev Quadrature
          DO i = 1, atom%n_shl
             ii = i
             IF (ang_loop == 1 ) THEN
                 ii=1
             END IF            
             write(iout,8) i
             Call Atomic_Grid(atom,atom%shl(ii)%ang%leb_ang,atom%shl(i),atom%rad%reg(i),center)
          END DO
!         Separable Quadrature
      ELSE
          DO i = 1, atom%n_shl
             ii = i
             IF (ang_loop == 1 ) THEN
                 ii=1
             END IF            
             atom%n_pts = atom%n_pts + atom%shl(ii)%ang%leb_ang%n_pts * atom%rad%reg(i)%n_pts
          END DO
          DO i = 1, atom%n_shl
             ii = i
             IF (ang_loop == 1 ) THEN
                 ii=1
             END IF       
             write(iout,9) i    
             Call Atomic_Grid(atom,atom%shl(ii)%ang%theta_ang,atom%shl(ii)%ang%theta_ang%reg,             &
                              atom%shl(ii)%ang%phi_ang,atom%shl(ii)%ang%phi_ang%reg,atom%rad%reg(i),      &
                              atom%shl(i),center)
          END DO
      END IF
      exact_vol = 4.d0 * pi * atom%R_max * atom%R_max * atom%R_max / 3.d0
      write(iout,*) 'Numerical Volume of Sphere = ', atom%vol
      write(iout,*) 'Exact Volume of Sphere     = ', exact_vol
      atom%local_time(3) = secnds(0.0) - t_start
      write(iout,*) 'Time for 3D grid generation = ',atom%local_time(3)
      atom%total_time(1) = atom%total_time(1) + atom%local_time(3)
  ELSE
      call lnkerr('need both radial and angular grids present to create 3D grid')
      Call exit      
  END IF
1  FORMAT(/,10x,' angular quadrature is non-separable (lebedev) in angles')
2  FORMAT(/,10x,'Lebedev quadrature data for shell = ',i4)
3  FORMAT(/,10x,'Radial quadrature data for shell = ',i4)
4  FORMAT(/,10x,' angular quadrature is separable in angles')
5  FORMAT(/,10x,'Theta quadrature data for shell = ',i4)
6  FORMAT(/,10x,'Phi quadrature data for shell = ',i4)
7  FORMAT(/,10x,'End Angular Coordinates for Shell = ',i4)
8  FORMAT(/,10x,'Computing 3D grid for shell = ',i2,' using Lebedev Angular Points')
9  FORMAT(/,10x,'Computing 3D grid for shell = ',i2,' using Gauss Angular Points')
END SUBROUTINE satshl
!********************************************************************************
!********************************************************************************
!deck Lebedev_Points
!***begin prologue     Lebedev_Points
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Theta grid
!***author             schneider, barry (nsf)
!***source             
!***purpose            Set boundary conditionsd for first and last shell
!***  
!***routines called
!***end prologue       Lebedev_Points
  SUBROUTINE Lebedev_Points(atom,leb,shl,rad_shl)
  IMPLICIT NONE
  TYPE(ATOMS)                                        :: atom
  TYPE(LEBEDEV)                                      :: leb
  TYPE(SHELLS)                                       :: shl
  TYPE(REGIONAL)                                     :: rad_shl
  shl%n_3d = rad_shl%n_pts * leb%n_pts
  atom%n_pts = atom%n_pts + shl%n_3d
  END SUBROUTINE Lebedev_Points
!********************************************************************************
!********************************************************************************
!deck Gauss_Points
!***begin prologue     Lebedev_Points
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Theta grid
!***author             schneider, barry (nsf)
!***source             
!***purpose            Set boundary conditionsd for first and last shell
!***  
!***routines called
!***end prologue       Gauss_Points
  SUBROUTINE Gauss_Points(atom,theta_ang,theta_reg,phi_ang,phi_reg,rad_shl,shl)
  IMPLICIT NONE
  TYPE(ATOMS)                                        :: atom
  TYPE(THETA)                                        :: theta_ang
  TYPE(PHI)                                          :: phi_ang
  TYPE(REGIONAL), DIMENSION(:)                       :: theta_reg
  TYPE(REGIONAL), DIMENSION(:)                       :: phi_reg
  TYPE(REGIONAL)                                     :: rad_shl
  TYPE(SHELLS)                                       :: shl
  INTEGER                                            :: i
  INTEGER                                            :: j
  shl%n_3d = 0
  DO i= 1, theta_ang%n_reg
     DO j = 1, phi_ang%n_reg
        shl%n_3d = shl%n_3d + theta_reg(i)%n_pts * phi_reg(j)%n_pts
     END DO
  END DO
  shl%n_2d=shl%n_3d
  shl%n_3d = shl%n_2d * rad_shl%n_pts
  atom%n_pts = atom%n_pts + shl%n_3d
  END SUBROUTINE Gauss_Points
!********************************************************************************
!********************************************************************************
!deck Set_BC
!***begin prologue     Set_BC
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Theta grid
!***author             schneider, barry (nsf)
!***source             
!***purpose            Set boundary conditionsd for first and last shell
!***  
!***routines called
!***end prologue       Set_BC
  SUBROUTINE Set_BC(coord_label,bc)
  IMPLICIT NONE
  TYPE(BOUNDARY_CONDITIONS)                 :: bc
  CHARACTER (*)                             :: coord_label
  INTEGER                                   :: i
  INTEGER                                   :: j
  INTEGER                                   :: nblock
  INTEGER                                   :: len
  CHARACTER (LEN=3)                         :: itoc
  INTEGER                                   :: intkey
  LOGICAL                                   :: logkey
  LOGICAL                                   :: dollar
  REAL(idp)                                 :: fpkey
!
  bc%n_fixed=intkey(card,'number_of_fixed_points',0,' ')
  bc%fixed(1)=.false.
  bc%fixed(2)=.false.
  bc%drop(1)=.false.
  bc%drop(2)=.false.
  IF(bc%n_fixed /= 0) THEN
     bc%fixed(1)=logkey(card,'left_fixed_point',.false.,' ')
     bc%fixed(2)=logkey(card,'right_fixed_point',.false.,' ')
     bc%drop(1)=logkey(card,'drop_left_function',.false.,' ')
     bc%drop(2)=logkey(card,'drop_right_function',.false.,' ')
  END IF
  write(iout,1) coord_label
  write(iout,2) bc%n_fixed,bc%fixed(1),bc%fixed(2),bc%drop(1),bc%drop(2)
1 Format(/,10x,'Coordinate Label = ',a16)
2 Format(/,15x,'number of fixed points    = ',i1,/15x,                   &
               'left point fixed          = ',l1,/,15x,                  &
               'right point fixed         = ',l1,/,15x,                  &
               'left point dropped        = ',l1,/15x,                   &
               'right point dropped       = ',l1)

  END SUBROUTINE Set_BC
!********************************************************************************
!********************************************************************************
!deck Fix_Angular_BC
!***begin prologue     Fix_Angular_BC
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           grid
!***author             schneider, barry (nsf)
!***source             
!***purpose            fix boundary conditions for first and last rgion
!***  
!***routines called
!***end prologue       Fix_BC
  SUBROUTINE Fix_Angular_BC(ang,reg,bc)
  IMPLICIT NONE
  TYPE(ANGULAR)                             :: ang
  TYPE(BOUNDARY_CONDITIONS)                 :: bc
  TYPE(REGIONAL), DIMENSION(:), ALLOCATABLE :: reg
  INTEGER                                   :: i
!
!                 Enforce boundary conditions.  They are taken to be independent of the specific
!                 shell but do need to be enforced on a shell level.
!
! Lets fix up the first region
!
! Special case of one region
!
  ang%n_phy=0
  DO i= 1, ang%n_reg
     ang%n_phy = ang%n_phy + reg(i)%n_pts - 2
  END DO
  ang%n_phy = ang%n_phy + 1
  IF (bc%drop(1) == .true.) THEN
      ang%n_phy = ang%n_phy - 1
  END IF
  IF (bc%drop(2) == .true.) THEN
      ang%n_phy = ang%n_phy - 1
  END IF
  IF (ang%n_reg == 1 ) THEN
      i = 1   
      IF (bc%drop(1) == .true.) THEN
          reg(i)%first = reg(i)%first + 1
          reg(i)%n_fun = reg(i)%n_fun - 1
      END IF
      IF (bc%drop(2) == .true.) THEN
          reg(i)%last = reg(i)%last - 1
          reg(i)%n_fun = reg(i)%n_fun - 1
      END IF
      IF (bc%n_fixed == 0 ) THEN
         reg(i)%type_quadrature = 'gauss'
         reg(i)%n_fixed = 0
      ELSE IF (bc%n_fixed == 1 ) THEN
         reg(i)%n_fixed = 1
         IF ( bc%fixed(1) == .true. ) THEN
              reg(i)%type_quadrature = 'radau'
         ELSE IF( bc%fixed(2) == .true. ) THEN
              reg(i)%type_quadrature = 'radau'
         END IF
      ELSE IF (bc%n_fixed == 2) THEN
         reg(i)%type_quadrature = 'lobatto'
         reg(i)%n_fixed = 2
      ELSE
         Call lnkerr('error in quadrature')
      END IF
  ELSE
      i = 1
      IF (bc%drop(1) == .true.) THEN
          reg(i)%first = reg(i)%first + 1
          reg(i)%n_fun = reg(i)%n_fun - 1
      END IF
      IF (bc%fixed(1) == .true. ) THEN
         reg(i)%type_quadrature = 'lobatto'
         reg(i)%n_fixed = 2
      ELSE
         reg(i)%type_quadrature = 'radau'
         reg(i)%n_fixed = 2
      END IF
      i = ang%n_reg
      IF (bc%drop(2) == .true.) THEN
          reg(i)%last = reg(i)%last - 1
          reg(i)%n_fun = reg(i)%n_fun - 1
      END IF
      IF (bc%fixed(2) == .true. )THEN
          reg(i)%type_quadrature = 'lobatto'
          reg(i)%n_fixed = 2
      ELSE
          reg(i)%type_quadrature = 'radau'
          reg(i)%n_fixed = 1
      END IF
  END IF
  END SUBROUTINE Fix_Angular_BC
!********************************************************************************
!********************************************************************************
!deck Fix_Radial_BC
!***begin prologue     Fix_Radial_BC
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           grid
!***author             schneider, barry (nsf)
!***source             
!***purpose            fix boundary conditions for first and last shell
!***  
!***routines called
!***end prologue       Fix_RAdial_BC
  SUBROUTINE Fix_Radial_BC(atom,shl,bc)
  IMPLICIT NONE
  TYPE(ATOMS)                               :: atom
  TYPE(BOUNDARY_CONDITIONS)                 :: bc
  TYPE(REGIONAL), DIMENSION(:), ALLOCATABLE :: shl
  INTEGER                                   :: i
!
!                 Enforce boundary conditions.  They are taken to be independent of the specific
!                 shell but do need to be enforced on a shell level.
!
! Get the total number of points with BC enforced
!
  atom%rad%n_phy=0
  DO i= 1, atom%n_shl
     atom%rad%n_phy = atom%rad%n_phy + shl(i)%n_pts - 2
  END DO
  atom%rad%n_phy = atom%rad%n_phy + 1
  IF (bc%drop(1) == .true.) THEN
      atom%rad%n_phy = atom%rad%n_phy - 1
  END IF
  IF (bc%drop(2) == .true.) THEN
      atom%rad%n_phy = atom%rad%n_phy - 1
  END IF
!
! Lets fix up the first region
!
! Special case of one region
!
  IF (atom%n_shl == 1 ) THEN
      i = 1   
      IF (bc%drop(1) == .true.) THEN
          shl(i)%first = shl(i)%first + 1
          shl(i)%n_fun = shl(i)%n_fun - 1
      END IF
      IF (bc%drop(2) == .true.) THEN
          shl(i)%last = shl(i)%last - 1
          shl(i)%n_fun = shl(i)%n_fun - 1
      END IF
      IF (bc%n_fixed == 0 ) THEN
         shl(i)%type_quadrature = 'gauss'
         shl(i)%n_fixed = 0
      ELSE IF (bc%n_fixed == 1 ) THEN
         shl(i)%n_fixed = 1
         IF ( bc%fixed(1) == .true. ) THEN
              shl(i)%type_quadrature = 'radau'
         ELSE IF( bc%fixed(2) == .true. ) THEN
              shl(i)%type_quadrature = 'radau'
         END IF
      ELSE IF (bc%n_fixed == 2) THEN
         shl(i)%type_quadrature = 'lobatto'
         shl(i)%n_fixed = 2
      ELSE
         Call lnkerr('error in quadrature')
      END IF
  ELSE
      i = 1
      IF (bc%drop(1) == .true.) THEN
          shl(i)%first = shl(i)%first + 1
          shl(i)%n_fun = shl(i)%n_fun - 1
      END IF
      IF (bc%fixed(1) == .true. ) THEN
         shl(i)%type_quadrature = 'lobatto'
         shl(i)%n_fixed = 2
      ELSE
         shl(i)%type_quadrature = 'radau'
         shl(i)%n_fixed = 1
      END IF
      i = atom%n_shl
      IF (bc%drop(2) == .true.) THEN
          shl(i)%last = shl(i)%last - 1
          shl(i)%n_fun = shl(i)%n_fun - 1
      END IF
      IF (bc%fixed(2) == .true. )THEN
          shl(i)%type_quadrature = 'lobatto'
          shl(i)%n_fixed = 2
      ELSE
          shl(i)%type_quadrature = 'radau'
          shl(i)%n_fixed = 1
      END IF
  END IF
  END SUBROUTINE Fix_Radial_BC
!********************************************************************************
!********************************************************************************
!deck Set_Angular_Variables
!***begin prologue     Set_Angular_Variables
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Angular grid
!***author             schneider, barry (nsf)
!***source             
!***purpose            Set boundary conditionsd for first and last shell
!***  
!***routines called
!***end prologue       Set_Angular_Variables
  SUBROUTINE Set_Angular_Variables(atom,ang,reg)
  IMPLICIT NONE
  TYPE(ANGULAR)                             :: ang
  TYPE(REGIONAL), DIMENSION(:), ALLOCATABLE :: reg
  TYPE(ATOMS)                               :: atom
  INTEGER                                   :: i
  INTEGER                                   :: j
  INTEGER                                   :: nblock
  INTEGER                                   :: len
  INTEGER                                   :: n_alloc
  CHARACTER (LEN=3)                         :: itoc
  CHARACTER (LEN=80)                        :: chrkey
  INTEGER                                   :: intkey
  LOGICAL                                   :: logkey
  LOGICAL                                   :: dollar
  REAL(idp)                                 :: fpkey
!
! Set the defaults for all regions.  The defaults are that lobatto quadrature is used, and both the left and right
! points of each element are fixed.  No basis functions are discarded. By enabling lobatto the Gauss routine does
! not even look at the variable n_fixed and simply knows that both endpoints are fiexed
! 
!
  automte=logkey(card,'automate',.false.,' ')
  reuse_sector_information=logkey(card,'reuse_sector_information',.false.,' ')
  IF(.not.automte) THEN
     ang%n_reg=intkey(card,'no_regions',1,' ')
     ALLOCATE(reg(1:ang%n_reg))
     DO i = 1, ang%n_reg
        CALL fparr(card,'region_'//itoc(i)//'_boundaries',reg(i)%edge,2,' ')
        reg(i)%n_pts = intkey(card,'region_'//itoc(i)//'_order',5,' ')
        reg(i)%type_quadrature = 'lobatto'
        reg(i)%n_fixed = 0 
        reg(i)%first = 1
        reg(i)%last  = reg(i)%n_pts
        reg(i)%n_fun = reg(i)%n_pts
    END DO
  ELSE
    defaults=logkey(card,'defaults',.false.,' ')
    IF (defaults == .true. ) THEN
        i=1
        nblock=1
        ang%n_reg=1
        nsubr=1
        nblock=1
        n_alloc=1
        boundl=-1.d0
        boundr=1.d0     
        IF(atom%coord_label(1:3) == 'phi') THEN
           boundl=0.d0
           boundr=two_pi   
           ang%phi_ang%fourier=logkey(card,'fourier',.false.,' ')
        END IF
        ALLOCATE(wa(1:n_alloc))
        call pakstr(str,len)
        IF ( dollar(str(1:len)//'_block_'//itoc(i),card, cpass,inp) ) THEN
             nply=intkey(card,'order',10,' ')
             wa(ang%n_reg)%edge(1) = boundl
             wa(ang%n_reg)%edge(2) = boundr
             wa(ang%n_reg)%n_pts = nply 
             write(iout,1) atom%coord_label, nply, boundl, boundr
        END IF
    ELSE
        ang%n_reg=0
        write(iout,*) '     Coordinate Label = ',atom%coord_label
        nblock=intkey(card,'number_of_major_blocks',1,' ')       
        ALLOCATE(wa(1:maxreg))
        call pakstr(str,len)
        DO i=1,nblock 
           IF ( dollar(str(1:len)//'_block_'//itoc(i),card, cpass,inp) ) THEN
                nsubr=intkey(card,'number_of_subregions',1,' ')
                skip=logkey(card,'skip',.false.,' ')
                IF(skip) THEN
                   write(iout,2) i
                ELSE 
                   nply=intkey(card,'order',10,' ')
                   nsubr=intkey(card,'number_of_subregions',1,' ')
                   boundl=fpkey(card,'left_boundary',0.d0,' ')
                   boundr=fpkey(card,'right_boundary',two_pi,' ')
                   step = ( boundr - boundl ) / nsubr
                   write(iout,3) i, nply, nsubr, boundl, boundr, step
                   DO j=1,nsubr
                      ang%n_reg= ang%n_reg + 1
                      wa(ang%n_reg)%edge(1) = boundl
                      wa(ang%n_reg)%edge(2) = wa(ang%n_reg)%edge(1) + step
                      boundl = wa(ang%n_reg)%edge(2)
                      wa(ang%n_reg)%n_pts = nply 
                   END DO
                END IF
           END IF
        END DO     
    END IF
    ALLOCATE(reg(1:ang%n_reg))
    ang%n_pts=0
    DO i = 1, ang%n_reg
       reg(i)%n_pts = wa(i)%n_pts
       ang%n_pts = ang%n_pts + reg(i)%n_pts
       reg(i)%edge(1:2) = wa(i)%edge(1:2)
       reg(i)%type_quadrature = 'lobatto'
       IF (ang%phi_ang%fourier == .true. ) THEN
           reg(i)%type_quadrature = 'trapezoidal'
       END IF
       reg(i)%n_fixed = 2 
       reg(i)%first = 1
       reg(i)%last  = reg(i)%n_pts
       reg(i)%n_fun = reg(i)%n_pts
    END DO
  END IF
  DEALLOCATE(wa)
1 FORMAT(/,1x, 'Using Default Values for = ',a8,         &
         /,15x,'quadrature order              = ',i4,    &
         /,15x,'left hand boundary            = ',e15.8, &
         /,15x,'right hand boundary           = ',e15.8)
2 FORMAT(/,1x,'skipping input block            = ',i4)
3 FORMAT(/,1x, 'block = ',i3,                            &
         /,15x,'quadrature order              = ',i4,    &
         /,15x,'number of subregions          = ',i4,    & 
         /,15x,'left hand boundary            = ',e15.8, &
         /,15x,'right hand boundary           = ',e15.8, &
         /,15x,'step                          = ',e15.8)
  END SUBROUTINE Set_Angular_Variables
!********************************************************************************
!********************************************************************************
!deck Set_Radial_Variables
!***begin prologue     Set_Radial_Variables
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Radial grid
!***author             schneider, barry (nsf)
!***source             
!***purpose            Set boundary conditionsd for first and last shell
!***  
!***routines called
!***end prologue       Set_Angular_Variables
  SUBROUTINE Set_Radial_Variables(atom,rad,shl)
  IMPLICIT NONE
  TYPE(ATOMS)                               :: atom
  TYPE(RADIAL)                              :: rad
  TYPE(REGIONAL), DIMENSION(:), ALLOCATABLE :: shl
  CHARACTER (LEN=80)                        :: chrkey
  INTEGER                                   :: intkey
  LOGICAL                                   :: logkey
  REAL(idp)                                 :: fpkey
  INTEGER                                   :: lenth 
  INTEGER                                   :: i
  INTEGER                                   :: len
  CHARACTER (LEN=3)                         :: itoc
  LOGICAL                                   :: dollar
!
! Set the defaults for all shells.  The defaults are that lobatto quadrature is used, and both the left and right
! points of each element are fixed.  No basis functions are discarded. By enabling lobatto the Gauss routine does
! not even look at the variable n_fixed and simply knows that both endpoints are fiexed
  ALLOCATE(shl(1:atom%n_shl))
  DO i = 1, atom%n_shl      ! Loop over the shells
     write(iout,1) i
     call pakstr(atom%coord_label,len)     
     str="$"//atom%coord_label(1:len)//'_shell_'//itoc(i)
     call pakstr(str,len)     
     IF ( dollar(str(1:len),card,cpass,inp) == .false. ) THEN  
          Call lnkerr('bad coordinate card')
     END IF
     defaults=logkey(card,'defaults',.false.,' ')!
     IF (defaults == .true. ) THEN
         shl(i)%n_pts = 10
         shl(i)%type_quadrature = chrkey(card,'type_quadrature','lobatto',' ' )
         shl(i)%edge(1)=0.d0
         shl(i)%edge(2)=10.d0
     ELSE
         shl(i)%n_pts=intkey(card,'order',10,' ')
         shl(i)%type_quadrature = chrkey(card,'type_quadrature','lobatto',' ' )
         shl(i)%edge(1)=fpkey(card,'left_boundary',0.d0,' ')
         shl(i)%edge(2)=fpkey(card,'right_boundary',10.d0,' ')
     END IF
     shl(i)%n_fun = shl(i)%n_pts
     shl(i)%first = 1
     shl(i)%last = shl(i)%n_pts
     shl(i)%n_fixed = 2 
     write(iout,2) shl(i)%n_pts, shl(i)%edge
  END DO
1  FORMAT(/,10x,'Radial quadrature data for shell = ',i4)
2 FORMAT(/,15x,'quadrature order              = ',i4,    &
         /,15x,'left hand boundary            = ',e15.8, &
         /,15x,'right hand boundary           = ',e15.8)
  END SUBROUTINE Set_Radial_Variables
!********************************************************************************
!********************************************************************************
!deck Lebedev_Grid.f
!***begin prologue     Lebedev_Grid
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Lebedev grid
!***author             schneider, barry (nsf)
!***source             
!***purpose            Points and weights of Lebedev quadrature.
!***  
!***routines called
!***end prologue       Lebedev_Grid
  SUBROUTINE Lebedev_Grid(ang)
  IMPLICIT NONE
  TYPE(LEBEDEV)                             :: ang
  INTEGER                                   :: rule_no
  INTEGER                                   :: intkey
  INTEGER                                   :: i
  LOGICAL                                   :: avail
!
  rule_no = intkey(card,'lebedev_rule_number',8,' ')
  ang%nang=rule_order_table(rule_no)
  avail=.false.
  IF ( rule_logic_table(rule_no) == 1 ) THEN
       avail = .true.
       precision = precision_table(rule_no)
       ang%nang=rule_order_table(rule_no)
       ang%nleb=precision
  ELSE
       Call lnkerr('error in lebedev rule')
  END IF
  WRITE(iout,1) ang%nleb, ang%nang
  CALL iosys ('write integer "lebedev quadrature order '//str//'" to grid',1,ang%nleb,0,' ')
  CALL iosys ('write integer "number lebedev points '//str//'" to grid',1,ang%nang,0,' ')
  atom%maxang=MAX(atom%maxang,ang%nang)
!
  order=ang%nang
  Call  Generate_Lebedev_Points_Weights(ang, leb_rule)
  Call Test_Rule(ang,rule_no)
!
  ang%sumwt = 0.d0
  DO i=1, ang%nang
     ang%sumwt = ang%sumwt + ang%w(i)
  END DO
  write(iout,2) ang%sumwt
  ang%w(:) = four_pi * ang%w(:)
1 FORMAT(/,10x,'Generating Lebedev Quadrature. Lebedev precision = ',i5,2x,   &
         'Number of lebedev points = ',i5)
2 FORMAT(/,10x,'Sum of Lebedev Weights = ',d15.8)
  END SUBROUTINE Lebedev_Grid
!********************************************************************************
!********************************************************************************
!deck Angular_Grid.f
!***begin prologue     Angular_Grid
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Angular grid
!***author             schneider, barry (nsf)
!***source             
!***purpose            read in angular information for theta based grid
!***  
!***routines called
!***end prologue       Angular_Grid
  SUBROUTINE Angular_Grid(atom,ang,reg)
  IMPLICIT NONE
  TYPE(ATOMS)                               :: atom
  TYPE(ANGULAR)                             :: ang
  TYPE(REGIONAL), DIMENSION(:), ALLOCATABLE :: reg
  TYPE(REAL_VECTOR)                         :: type_real_vector
  TYPE(REAL_MATRIX)                         :: type_real_matrix
  INTEGER                                   :: intkey
  INTEGER                                   :: i
  INTEGER                                   :: j
  CHARACTER (LEN=30)                        :: chrkey
  CHARACTER (LEN=3)                         :: itoc
!
! We define a single theta grid independent of the m quantum number.  The grid needs to capture the
! polynomic behavior of all P(l,m) on the interval {-1,1} with a factor of ( 1-x*x)^{1/2} removed.  With the factor 
! removed, these functions are all polynomials in {x}.  How well the exact behavior is captured depends on the 
! size of the quadrature.
!
  DO i= 1, ang%n_reg 
     Write(iout,1) i, reg(i)%n_pts, reg(i)%type_quadrature, reg(i)%n_fixed
     ALLOCATE( reg(i)%q(1:reg(i)%n_pts), reg(i)%wt(1:reg(i)%n_pts) )
     Call gauss(reg(i)%q, reg(i)%wt, reg(i)%edge,type_quadrature=reg(i)%type_quadrature,                     &
                fixed_point=reg(i)%n_fixed, n=reg(i)%n_pts)
  END DO
  IF (iosys_on == .true.) THEN
      DO i= 1, ang%n_reg 
         CALL iosys('write integer "no '//atom%coord_label//' quadrature points '//str//' region '               &
                    //itoc(i)//'" to grid',1,reg(i)%n_pts,0,' ')
         CALL iosys('write character "'//atom%coord_label//' quadrature type '//str//' region '                  &
                    //itoc(i)//'"to grid',0,0,0,reg(i)%type_quadrature)
         CALL iosys('write real "'//atom%coord_label//' quadrature points '//str//' region '                     &
                     //itoc(i)//'" to grid',reg(i)%n_pts,reg(i)%q,0,' ')
         CALL iosys('write real "'//atom%coord_label//' quadrature weights '//str//' region '                    &
                    //itoc(i)//'" to grid',reg(i)%n_pts,reg(i)%wt,0,' ')
      END DO
  END IF
!
  IF (print_ang(1) == .true.) THEN
      DO i = 1, ang%n_reg
         call Print_Matrix(type_real_vector,reg(i)%q,title='Points Region-'//itoc(i))
         call Print_Matrix(type_real_vector,reg(i)%wt,title='Weights Region-'//itoc(i))
      END DO
  END IF
!
  IF( fedvr_functions == .true. ) THEN
      write(iout,*) ' Computing FEDVR functions'
!     Calculate factors needed for matrix elements
!
      IF (atom%coord_label == 'theta' ) THEN
          DO i = 1, ang%n_reg
             ALLOCATE( reg(i)%q_fac(1:reg(i)%n_pts), reg(i)%inv_q_fac(1:reg(i)%n_pts),                               &
                       reg(i)%inv_sqrt_q_fac(1:reg(i)%n_pts))
             reg(i)%q_fac(:) = one - reg(i)%q(:) * reg(i)%q(:)
             reg(i)%inv_q_fac(:) = one / reg(i)%q_fac(:)
             reg(i)%inv_sqrt_q_fac(:) = Sqrt ( reg(i)%inv_q_fac(:) )
          END DO
      ELSE
          DO i = 1, ang%n_reg
             ALLOCATE( reg(i)%q_fac(1:reg(i)%n_pts), reg(i)%inv_q_fac(1:reg(i)%n_pts),                               &
                       reg(i)%inv_sqrt_q_fac(1:reg(i)%n_pts))
             reg(i)%q_fac(:) = 1.d0
             reg(i)%inv_q_fac(:) = 1.d0
             reg(i)%inv_sqrt_q_fac(:) = 1.d0
          END DO
      END IF
!
      IF (print_ang(3) == .true.) THEN
          DO i = 1, ang%n_reg
             call Print_Matrix(type_real_vector,reg(i)%q_fac,title='q_factor Region-'//itoc(i))
             call Print_Matrix(type_real_vector,reg(i)%inv_q_fac,title='inverse_q_factor Region-'//itoc(i))
             call Print_Matrix(type_real_vector,reg(i)%inv_sqrt_q_fac,title='inverse_sqrt_q_factor Region-'//itoc(i))
          END DO
!
      END IF  
!
      DO i= 1, ang%n_reg 
         ALLOCATE( reg(i)%p(1:reg(i)%n_pts,1:reg(i)%n_pts),reg(i)%dp(1:reg(i)%n_pts,1:reg(i)%n_pts),               &
                   reg(i)%ddp(1:reg(i)%n_pts,1:reg(i)%n_pts),reg(i)%normalization(1:reg(i)%n_pts) )
         CALL cpoly(reg(i)%p,reg(i)%dp,reg(i)%ddp,reg(i)%q,reg(i)%n_pts)
      END DO
!
      IF (print_ang(2) == .true.) THEN
          DO i = 1, ang%n_reg
             call Print_Matrix(type_real_matrix,reg(i)%p,reg(i)%n_pts,reg(i)%n_pts,                     &
                               title='Unnormalized Polynomials Region-'//itoc(i))
             call Print_Matrix(type_real_matrix,reg(i)%dp,reg(i)%n_pts,reg(i)%n_pts,                    &
                               title='First Derivative of Unnormalized Polynomials Region-'//itoc(i))
             call Print_Matrix(type_real_matrix,reg(i)%ddp,reg(i)%n_pts,reg(i)%n_pts,                   &
                               title='Second Derivative of Unnormalized Polynomials Region-'//itoc(i))
         END DO
      END IF
!
! Get normalization factors and renormalize
!
      Call Normalization(reg,ang%n_reg)
!                                                                                                                
!                  Normalize the functions in their own region.  
      DO i = 1, ang%n_reg
         ALLOCATE( reg(i)%bf(1:reg(i)%n_pts,1:reg(i)%n_pts), reg(i)%dbf(1:reg(i)%n_pts,1:reg(i)%n_pts), &
                   reg(i)%ddbf(1:reg(i)%n_pts,1:reg(i)%n_pts)) 
         DO j = 1, reg(i)%n_pts
            reg(i)%bf(:,j)   = reg(i)%p(:,j) / sqrt( reg(i)%wt(j) )
            reg(i)%dbf(:,j)  = reg(i)%dp(:,j) / sqrt( reg(i)%wt(j) )
            reg(i)%ddbf(:,j) = reg(i)%ddp(:,j) / sqrt( reg(i)%wt(j) )
         END DO
      END DO
!
      IF (print_ang(2) == .true.) THEN
          DO i = 1, ang%n_reg
             call Print_Matrix(type_real_matrix,reg(i)%bf(:,reg(i)%first:reg(i)%last),                  &
                               reg(i)%n_pts,reg(i)%n_fun,                                               &
                               title='Final Region Normalized Basis Polynomials Region-'//itoc(i))
             call Print_Matrix(type_real_matrix,reg(i)%dbf(:,reg(i)%first:reg(i)%last),                 &
                               reg(i)%n_pts,reg(i)%n_fun,                                               &
                               title='First Derivative of Final Region Normalized Basis Polynomials Region-'//itoc(i))
             call Print_Matrix(type_real_matrix,reg(i)%ddbf(:,reg(i)%first:reg(i)%last),                &
                               reg(i)%n_pts,reg(i)%n_fun,                                               &
                               title='Second Derivative of Final Region Normalized Basis Polynomials Region-'//itoc(i))
         END DO
      END IF  
      IF (iosys_on == .true.) THEN
          DO i= 1, ang%n_reg 
             CALL iosys('write real "'//atom%coord_label//' q_fac '//str//' region '//itoc(i)//'" to grid',          &
                         reg(i)%n_pts,reg(i)%q_fac,0,' ')
             CALL iosys('write real "'//atom%coord_label//' inv_q_fac '//str//' region '//itoc(i)//'" to grid',      &
                         reg(i)%n_pts,reg(i)%inv_q_fac,0,' ')
             CALL iosys('write real "'//atom%coord_label//' inv_sqrt_q_fac '//str//' region '//itoc(i)//'" to grid', &
                         reg(i)%n_pts,reg(i)%inv_sqrt_q_fac,0,' ')
             CALL iosys('write real "'//atom%coord_label//' polynomials '//str//' region '//itoc(i)//'" to grid',      &
                         reg(i)%n_pts,reg(i)%p,0,' ')
             CALL iosys('write real "first derivative of '//atom%coord_label//' polynomials '//str//' region '         &
                         //itoc(i)//'" to grid',reg(i)%n_pts,reg(i)%dp,0,' ')
             CALL iosys('write real "second derivative of '//atom%coord_label//' polynomials '//str//' region '        &
                        //itoc(i)//'" to grid',reg(i)%n_pts,reg(i)%ddp,0,' ')
             CALL iosys('write real "'//atom%coord_label//' normalized polynomials '//str//' region '       &
                        //itoc(i)//'" to grid',reg(i)%n_pts,reg(i)%bf,0,' ')
             CALL iosys('write real "first derivative of normalized '//atom%coord_label//' polynomials '    &
                        //str//' region '//itoc(i)//'" to grid',reg(i)%n_pts,reg(i)%dbf,0,' ')
             CALL iosys('write real "second derivative of normalized '//atom%coord_label//' polynomials '   &
                        //str//' region '//itoc(i)//'" to grid',reg(i)%n_pts,reg(i)%ddbf,0,' ')
         END DO
      END IF
  ELSE
      write(iout,*) ' Not computing FEDVR functions'
  END IF
1 FORMAT(/,' Generating Angular Quadrature Region = ',i3,2x,                         &
           ' Number of Points = ',i5,/,' Type Quadrature = ',a8,2x,'left(1) or right(2) fixed point = ',i1)
  END SUBROUTINE Angular_Grid
!********************************************************************************
!********************************************************************************
!deck Radial_Grid.f
!***begin prologue     Radial_Grid
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Radial_Grid, link m6200
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            read in radial quadrature information.
!***routines called
!***end prologue       Radial_Grid
  SUBROUTINE Radial_Grid(atom,rad,shl)
  IMPLICIT NONE
  TYPE(ATOMS)                                  :: atom
  TYPE(RADIAL)                                 :: rad
  TYPE(REGIONAL), DIMENSION(:)                 :: shl
  TYPE(REAL_VECTOR)                            :: type_real_vector
  TYPE(REAL_MATRIX)                            :: type_real_matrix
  INTEGER                                      :: intkey
  INTEGER                                      :: i
  INTEGER                                      :: j
  INTEGER                                      :: len
  LOGICAL                                      :: dollar
  CHARACTER (LEN=30)                           :: chrkey
  CHARACTER (LEN=3)                            :: itoc
!
  DO i=1, atom%n_shl
     str="$"//atom%coord_label(1:len)//'_shell_'//itoc(i)
     Write(iout,1) i, shl(i)%n_pts, shl(i)%type_quadrature, shl(i)%n_fixed
     ALLOCATE( shl(i)%q(1:shl(i)%n_pts), shl(i)%wt(1:shl(i)%n_pts) )
     Call gauss(shl(i)%q, shl(i)%wt, shl(i)%edge,type_quadrature=shl(i)%type_quadrature,                     &
                fixed_point=shl(i)%n_fixed, n=shl(i)%n_pts)
  END DO
  IF (iosys_on == .true. ) THEN
      DO i=1, atom%n_shl
         CALL iosys('write integer "no '//atom%coord_label//' quadrature points '//str//'" to grid',             &
                     1,shl(i)%n_pts,0,' ')
         CALL iosys('write character "'//atom%coord_label//' quadrature type '//str//'"to grid',                 &
                     0,0,0,shl(i)%type_quadrature)
         CALL iosys('write real "'//atom%coord_label//' quadrature points '//str//'" to grid',                   &
                     shl(i)%n_pts,shl(i)%q,0,' ')
         CALL iosys('write real "'//atom%coord_label//' quadrature weights '//str//'" to grid',                  &
                     shl(i)%n_pts,shl(i)%wt,0,' ')
      END DO
  END IF
!
  IF (print_rad(1) == .true.) THEN
      DO i = 1, atom%n_shl
         call Print_Matrix(type_real_vector,shl(i)%q,title='Points Shell-'//itoc(i))
         call Print_Matrix(type_real_vector,shl(i)%wt,title='Weights Shell-'//itoc(i))
     END DO
  END IF
  IF( fedvr_functions == .true. ) THEN
      write(iout,*) ' Computing FEDVR functions'
!
!     Calculate factors needed for matrix elements
!
      DO i = 1, atom%n_shl  
         ALLOCATE( shl(i)%q_fac(1:shl(i)%n_pts), shl(i)%inv_q_fac(1:shl(i)%n_pts), shl(i)%inv_sqrt_q_fac(1:shl(i)%n_pts))
         shl(i)%q_fac(:) = 1.d0
         shl(i)%inv_q_fac(:) = 1.d0
         shl(i)%inv_sqrt_q_fac(:) = 1.d0
      END DO
!
      IF (print_rad(3) == .true.) THEN
          DO i = 1, atom%n_shl  
             call Print_Matrix(type_real_vector,shl(i)%q_fac,title='q_factor Shell-'//itoc(i))
             call Print_Matrix(type_real_vector,shl(i)%inv_q_fac,title='inverse_q_factor Shell-'//itoc(i))
             call Print_Matrix(type_real_vector,shl(i)%inv_sqrt_q_fac,title='inverse_sqrt_q_factor Shell-'//itoc(i))
         END DO  
      END IF
!
      DO i = 1, atom%n_shl  
         ALLOCATE( shl(i)%p(1:shl(i)%n_pts,1:shl(i)%n_pts),shl(i)%dp(1:shl(i)%n_pts,1:shl(i)%n_pts),                &
                   shl(i)%ddp(1:shl(i)%n_pts,1:shl(i)%n_pts),shl(i)%normalization(1:shl(i)%n_pts) )
         CALL cpoly(shl(i)%p,shl(i)%dp,shl(i)%ddp,shl(i)%q,shl(i)%n_pts)
      END DO
      IF (print_rad(2) == .true.) THEN
          DO i = 1, atom%n_shl  
             call Print_Matrix(type_real_matrix,shl(i)%p,shl(i)%n_pts,shl(i)%n_pts,                                 &
                               title='Unnormalized Polynomials Shell-'//itoc(i))
             call Print_Matrix(type_real_matrix,shl(i)%dp,shl(i)%n_pts,shl(i)%n_pts,                                &
                               title='First Derivative of Unnormalized Polynomials Shell-'//itoc(i))
             call Print_Matrix(type_real_matrix,shl(i)%ddp,shl(i)%n_pts,shl(i)%n_pts,                               &
                               title='Second Derivative of Unnormalized Polynomials Shell-'//itoc(i))
          END DO
      END IF
!
!     Get normalization factors and renormalize
!
      Call Normalization(shl,atom%n_shl)
!                                                                                                                
!                  Normalize the functions in their own region.  
      DO i = 1, atom%n_shl  
         ALLOCATE( shl(i)%bf(1:shl(i)%n_pts,1:shl(i)%n_pts), shl(i)%dbf(1:shl(i)%n_pts,1:shl(i)%n_pts), &
                   shl(i)%ddbf(1:shl(i)%n_pts,1:shl(i)%n_pts)) 
         DO j = 1, shl(i)%n_pts
            shl(i)%bf(:,j)   = shl(i)%p(:,j) / sqrt( shl(i)%wt(j) )
            shl(i)%dbf(:,j)  = shl(i)%dp(:,j) / sqrt( shl(i)%wt(j) )
            shl(i)%ddbf(:,j) = shl(i)%ddp(:,j) / sqrt( shl(i)%wt(j) )
         END DO
      END DO
!
      IF (print_rad(2) == .true.) THEN
          DO i = 1, atom%n_shl   
             call Print_Matrix(type_real_matrix,shl(i)%bf(:,shl(i)%first:shl(i)%last),                  &
                               shl(i)%n_pts,shl(i)%n_fun,                                               &
                               title='Final Region Normalized Basis Polynomials Region-'//itoc(i))
             call Print_Matrix(type_real_matrix,shl(i)%dbf(:,shl(i)%first:shl(i)%last),                 &
                               shl(i)%n_pts,shl(i)%n_fun,                                               &
                               title='First Derivative of Final Region Normalized Basis Polynomials Region-'//itoc(i))
             call Print_Matrix(type_real_matrix,shl(i)%ddbf(:,shl(i)%first:shl(i)%last),                &
                               shl(i)%n_pts,shl(i)%n_fun,                                               &
                               title='Second Derivative of Final Region Normalized Basis Polynomials Region-'//itoc(i))
         END DO
      END IF
      IF (iosys_on == .true.) THEN
          DO i=1, atom%n_shl
             CALL iosys('write real "'//atom%coord_label//' q_fac '//str//'" to grid',          &
                         shl(i)%n_pts,shl(i)%q_fac,0,' ')
             CALL iosys('write real "'//atom%coord_label//' inv_q_fac '//str//'" to grid',      &
                         shl(i)%n_pts,shl(i)%inv_q_fac,0,' ')
             CALL iosys('write real "'//atom%coord_label//' inv_sqrt_q_fac '//str//'" to grid', &
                         shl(i)%n_pts,shl(i)%inv_sqrt_q_fac,0,' ')
             CALL iosys('write real "'//atom%coord_label//' polynomials '//str//'" to grid',                            &
                         shl(i)%n_pts,shl(i)%p,0,' ')
             CALL iosys('write real "first derivative of '//atom%coord_label//' polynomials '//str//'" to grid',        &
                         shl(i)%n_pts,shl(i)%dp,0,' ')
             CALL iosys('write real "second derivative of '//atom%coord_label//' polynomials '//str//'" to grid',       &
                         shl(i)%n_pts,shl(i)%ddp,0,' ')
             CALL iosys('write real "'//atom%coord_label//' normalized polynomials '//str//' region '       &
                        //itoc(i)//'" to grid',shl(i)%n_pts,shl(i)%bf,0,' ')
             CALL iosys('write real "first derivative of normalized '//atom%coord_label//' polynomials '    &
                        //str//' region '//itoc(i)//'" to grid',shl(i)%n_pts,shl(i)%dbf,0,' ')
             CALL iosys('write real "second derivative of normalized '//atom%coord_label//' polynomials '   &
                        //str//' region '//itoc(i)//'" to grid',shl(i)%n_pts,shl(i)%ddbf,0,' ')
         END DO
      END IF
  ELSE
      write(iout,*) ' Not computing FEDVR functions'
  END IF  
!
1 FORMAT(' radial region = ',i3,2x,' number of points = ',i5,2x,' type quadrature = ',a8,           &
         ' number of fixed points = ',i1)
  END SUBROUTINE Radial_Grid
!***********************************************************************
!***********************************************************************
!deck Three_d_Lebedev_Grid
!***begin prologue    Three_d_Lebedev_Grid
!***date written       000702   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr
!***
!***author             schneider, b. i.(nsf)
!***source             
!***purpose            Compute the cartesian and spherical grids.  The input are the
!***                   points and weights defined with respect to each atomic center.
!***                   The subroutine handles only a single shel not the entire center.
!***                   The output are those same quantities define with respect to a fixed
!***                   center taken at the origin.
!***description        
!***                   
!***                   
!***                   
!***                   

!***references         

!***routines called    iosys, util and mdutil
!***end prologue       Cartesion_Grid
  Subroutine Three_d_Lebedev_Grid(atom,leb,shl,rad_shl,center)
  IMPLICIT NONE
  TYPE(ATOMS)                                        :: atom
  TYPE(LEBEDEV)                                      :: leb
  TYPE(SHELLS)                                       :: shl
  TYPE(ATOMS), DIMENSION(:)                          :: center
  TYPE(REGIONAL)                                     :: rad_shl
  TYPE(REAL_MATRIX)                                  :: type_real_matrix
  TYPE(REAL_VECTOR)                                  :: type_real_vector
  CHARACTER(LEN=1), DIMENSION(7)                     :: collab
  REAL(idp)                                          :: atan2
  REAL(idp)                                          :: sdot
  REAL(idp)                                          :: dist
  INTEGER                                            :: i
  INTEGER                                            :: nc
  INTEGER                                            :: start
  INTEGER                                            :: end
  shl%n_3d = rad_shl%n_pts * leb%n_pts
  atom%n_pts = atom%n_pts + shl%n_3d
  write(iout,1) shl%n_3d
  ALLOCATE(shl%xyz_grid(1:shl%n_3d,1:3), shl%weight(1:shl%n_3d), shl%sph_grid(1:shl%n_3d,1:8),  &
           shl%integration_weight(1:shl%n_3d) )
  start = 0
  end = 0
!
!     The coodinates were all computed with respect to the atomic center.  So, what we have
!     stored is x_i, y_i and z_i.  Now get the coordinates with respect to a fixed origin. 
!     That is why we ADD not SUBTRACT the atomic center positions.  The grid is now defined with
!     respect to a single origin.
  
  DO i = 1, rad_shl%n_pts
     start = end + 1
     end = end + leb%n_pts
     shl%weight(start:end)     = rad_shl%wt(i) * leb%w(:) ! total weight
     shl%integration_weight(start:end)    = rad_shl%wt(i) * rad_shl%q(i) * rad_shl%q(i) * leb%w(:) ! composite weight
     shl%xyz_grid(start:end,1) = rad_shl%q(i) * leb%q(1,:) + atom%cen(1) 
                                          ! x = x_i + X_i
     shl%xyz_grid(start:end,2) = rad_shl%q(i) * leb%q(2,:) + atom%cen(2)
                                          ! y = y_i + Y_i
     shl%xyz_grid(start:end,3) = rad_shl%q(i) * leb%q(3,:) + atom%cen(3)                 
                                          ! z = z_i + Z_i
     shl%sph_grid(start:end,8) = shl%xyz_grid(start:end,1) * shl%xyz_grid(start:end,1) + &
                                 shl%xyz_grid(start:end,2) * shl%xyz_grid(start:end,2) + &
                                 shl%xyz_grid(start:end,3) * shl%xyz_grid(start:end,3)
                                          ! r*r
     shl%sph_grid(start:end,1) = sqrt ( shl%xyz_grid(start:end,8) )   
                                          ! r
     shl%sph_grid(start:end,2) = shl%xyz_grid(start:end,3) / shl%sph_grid(start:end,1)
                                          ! z/r = cos(theta)
     shl%sph_grid(start:end,3) = acos ( shl%sph_grid(start:end,2) )
                                          ! theta
     shl%sph_grid(start:end,4) = sqrt ( 1.d0 - shl%sph_grid(start:end,2) * shl%sph_grid(start:end,2) )
                                          ! sin(theta)
!     shl%sph_grid(start:end,5) = 1.d0
!     shl%sph_grid(start:end,6) = 0.d0
     shl%sph_grid(start:end,7) = atan2( shl%xyz_grid(start:end,2),shl%xyz_grid(start:end,1) )     
     shl%sph_grid(start:end,5) = cos( shl%sph_grid(start:end,7) )     
!     shl%sph_grid(start:end,6) = sin( shl%sph_grid(start:end,7) )     
     shl%sph_grid(start:end,6) = sqrt ( 1.d0 - shl%sph_grid(start:end,5) * shl%sph_grid(start:end,5) )
!     shl%sph_grid(start:end,6) = sin( shl%sph_grid(start:end,7) )     
!     DO j = start, end
!        IF ( shl%sph_grid(j,4) /= 0.d0) THEN
!             shl%sph_grid(j,5) = shl%xyz_grid(j,1) / (shl%sph_grid(j,1) * shl%sph_grid(j,4) )
!                                          ! cos(phi)
!             shl%sph_grid(j,6) = shl%xyz_grid(j,2) / (shl%sph_grid(j,1) * shl%sph_grid(j,4) )
!                                          ! sin(phi)
!        END IF
!     END DO 
!    shl%sph_grid(start:end,7) = atan2( shl%sph_grid(start:end,6),shl%sph_grid(start:end,5) )
                                          ! phi
  END DO
  atom%vol = atom%vol + sdot(shl%n_3d,shl%weight(start:),1,shl%sph_grid(start:,8),1)
  IF (print_leb(1) == .true. ) THEN
      collab(1)='x'
      collab(2)='y'
      collab(3)='z'
      Call Print_Matrix(type_real_matrix,shl%xyz_grid,shl%n_3d,3,          &
                        title='Cartesian Points',frmt='fr',collab=collab)
      Call Print_Matrix(type_real_vector,shl%weight,                       &
                        frmt='fr',title='Cartesian Weights')
  END IF
  IF (print_leb(2) == .true. ) THEN
      Call Print_Matrix(type_real_vector,shl%sph_grid(:,1),title='r',frmt='fr')
      Call Print_Matrix(type_real_vector,shl%sph_grid(:,2),title='Cos theta Points',frmt='fr')
      Call Print_Matrix(type_real_vector,shl%sph_grid(:,3),title='Theta Points',frmt='fr')
      Call Print_Matrix(type_real_vector,shl%sph_grid(:,4),title='Sin Theta Points',frmt='fr')
      Call Print_Matrix(type_real_vector,shl%sph_grid(:,5),title='Cos Phi Points',frmt='fr')
      Call Print_Matrix(type_real_vector,shl%sph_grid(:,6),title='Sin Phi Points',frmt='fr')
      Call Print_Matrix(type_real_vector,shl%sph_grid(:,7),title='Phi Points',frmt='fr')
  END IF
  IF (yukawa_on == .true. ) THEN
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
  END IF
1 FORMAT(/,10x,'Number of 3D points for this shell = ',i8)
END SUBROUTINE Three_d_Lebedev_grid
!***********************************************************************
!***********************************************************************
!deck Three_d_Gauss_Grid
!***begin prologue     Three_d_Gauss_Grid
!***date written       000702   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr
!***
!***author             schneider, b. i.(nsf)
!***source             
!***purpose            
!***                   
!***description        
!***                   
!***                   
!***                   
!***                   

!***references         

!***routines called    iosys, util and mdutil
!***end prologue       Cartesion_Grid
  Subroutine Three_d_Gauss_Grid(atom,theta_ang,theta_reg,phi_ang,phi_reg,rad_shl,shl,center)
  IMPLICIT NONE
  TYPE(ATOMS)                                        :: atom
  TYPE(THETA)                                        :: theta_ang
  TYPE(PHI)                                          :: phi_ang
  TYPE(REGIONAL), DIMENSION(:)                       :: theta_reg
  TYPE(REGIONAL), DIMENSION(:)                       :: phi_reg
  TYPE(REGIONAL)                                     :: rad_shl
  TYPE(SHELLS)                                       :: shl
  TYPE(ATOMS), DIMENSION(:)                          :: center
  REAL(idp), DIMENSION(:,:), ALLOCATABLE             :: temp
  INTEGER                                            :: i
  INTEGER                                            :: j
  INTEGER                                            :: ii
  INTEGER                                            :: jj
  INTEGER                                            :: nc
  CHARACTER(LEN=1), DIMENSION(7)                     :: collab
  REAL(idp)                                          :: atan2
  REAL(idp)                                          :: sdot
  REAL(idp)                                          :: dist
  INTEGER                                            :: start
  INTEGER                                            :: end
  INTEGER                                            :: count
  INTEGER                                            :: n_2d
  DO i = 1, theta_ang%n_reg
     ALLOCATE(theta_reg(i)%sin_theta(1:theta_reg(i)%n_pts))
     theta_reg(i)%sin_theta(:) = sqrt ( 1.d0 - theta_reg(i)%q(:) * theta_reg(i)%q(:) )
  END DO
  DO i = 1, phi_ang%n_reg
     ALLOCATE(phi_reg(i)%sin_phi(1:phi_reg(i)%n_pts), phi_reg(i)%cos_phi(1:phi_reg(i)%n_pts))
     phi_reg(i)%cos_phi(:) = cos( phi_reg(i)%q(:) )     
     phi_reg(i)%sin_phi(:) = sqrt ( 1.d0 - phi_reg(i)%cos_phi(:) * phi_reg(i)%cos_phi(:) )
  END DO
  shl%n_3d = 0
  DO i= 1, theta_ang%n_reg
     DO j = 1, phi_ang%n_reg
        shl%n_3d = shl%n_3d + theta_reg(i)%n_pts * phi_reg(j)%n_pts
     END DO
  END DO
  shl%n_2d=shl%n_3d
  shl%n_3d = shl%n_2d * rad_shl%n_pts
  atom%n_pts = atom%n_pts + shl%n_3d
  write(iout,1) shl%n_3d
  ALLOCATE ( temp(1:n_2d,1:4) )
  count = 0
  DO i = 1, theta_ang%n_reg
     DO j =1, phi_ang%n_reg
        DO ii = 1, theta_reg(i)%n_pts
           DO jj = 1, phi_reg(j)%n_pts
              count = count + 1
              temp(count,1) = theta_reg(i)%wt(ii) * phi_reg(j)%wt(jj) ! total angular weight
              temp(count,2) = theta_reg(i)%sin_theta(ii) * phi_reg(j)%cos_phi(jj)
              temp(count,3) = theta_reg(i)%sin_theta(ii) * phi_reg(j)%sin_phi(jj)
              temp(count,4) = theta_reg(i)%q(ii) 
           END DO
        END DO
     END DO
  END DO
  ALLOCATE(shl%xyz_grid(1:shl%n_3d,1:3), shl%weight(1:shl%n_3d), shl%integration_weight(1:shl%n_3d))
  start = 0
  end = 0
  DO i = 1, rad_shl%n_pts 
     start = end + 1
     end = end + n_2d
     shl%weight(start:end) =     rad_shl%wt(i) * temp(:,1) ! total weight
     shl%integration_weight(start:end)    = rad_shl%wt(i) * rad_shl%q(i) * rad_shl%q(i) * temp(:,1) ! composite weight
     shl%xyz_grid(start:end,1) = rad_shl%q(i)  * temp(:,2) + atom%cen(1) 
                                          ! x = x_i + X_i
     shl%xyz_grid(start:end,2) = rad_shl%q(i)  * temp(:,3) + atom%cen(2)
                                          ! y = y_i + Y_i
     shl%xyz_grid(start:end,3) = rad_shl%q(i)  * temp(:,4) + atom%cen(3)                 
                                          ! z = z_i + Z_i
  END DO
  DEALLOCATE ( temp )
  ALLOCATE(shl%sph_grid(1:shl%n_3d,1:6))
  start = 0
  end = 0
  DO i = 1, rad_shl%n_pts 
     start = end + 1
     end = end + n_2d
     shl%sph_grid(start:end,1) = sqrt ( shl%xyz_grid(start:end,1) * shl%xyz_grid(start:end,1) + &
                                        shl%xyz_grid(start:end,2) * shl%xyz_grid(start:end,2) + &
                                        shl%xyz_grid(start:end,3) * shl%xyz_grid(start:end,3) )   
                                          ! r
     shl%sph_grid(start:end,2) = 0.d0
                                          ! z/r = cos(theta)
     DO j = start, end
        IF ( shl%sph_grid(j,1)  /= 0.d0 ) THEN
             shl%sph_grid(j,2) = shl%xyz_grid(j,3) / shl%sph_grid(j,1)
                                              ! z/r = cos(theta)
        END IF
     END DO
     shl%sph_grid(start:end,3) = acos ( shl%sph_grid(start:end,2) )
                                              ! theta
     shl%sph_grid(start:end,4) = sqrt ( 1.d0 - shl%sph_grid(start:end,2) * shl%sph_grid(start:end,2) )
                                              ! sin(theta)
     shl%sph_grid(start:end,7) = atan2( shl%xyz_grid(start:end,2),shl%xyz_grid(start:end,1) )     
     shl%sph_grid(start:end,5) = cos( shl%sph_grid(start:end,7) )     
     shl%sph_grid(start:end,6) = sqrt ( 1.d0 - shl%sph_grid(start:end,5) * shl%sph_grid(start:end,5) )
  END DO
  DO i = 1, shl%n_3d 
     atom%vol  = atom%vol  + shl%weight(i) * shl%sph_grid(i,1) * shl%sph_grid(i,1)
  END DO
  IF (print_gauss(1) == .true. ) THEN
      collab(1)='x'
      collab(2)='y'
      collab(3)='z'
      Call Print_Matrix(type_real_matrix,shl%xyz_grid,shl%n_3d,3,          &
                        title='Cartesian Points',frmt='fr',collab=collab)
      Call Print_Matrix(type_real_vector,shl%weight,                       &
                        frmt='fr',title='Cartesian Weights')
  END IF
  IF (print_gauss(2) == .true. ) THEN
      Call Print_Matrix(type_real_vector,shl%sph_grid(:,1),title='r',frmt='fr')
      Call Print_Matrix(type_real_vector,shl%sph_grid(:,2),title='Cos theta Points',frmt='fr')
      Call Print_Matrix(type_real_vector,shl%sph_grid(:,3),title='Theta Points',frmt='fr')
      Call Print_Matrix(type_real_vector,shl%sph_grid(:,4),title='Sin Theta Points',frmt='fr')
      Call Print_Matrix(type_real_vector,shl%sph_grid(:,5),title='Cos Phi Points',frmt='fr')
      Call Print_Matrix(type_real_vector,shl%sph_grid(:,6),title='Sin Phi Points',frmt='fr')
      Call Print_Matrix(type_real_vector,shl%sph_grid(:,7),title='Phi Points',frmt='fr')
  END IF
  IF (yukawa_on == .true. ) THEN
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
  END IF
1 FORMAT(/,10x,'Number of 3D points for this shell = ',i8)
END SUBROUTINE Three_d_Gauss_Grid
!***********************************************************************
!***********************************************************************
END MODULE Grid_Generation
!***********************************************************************
!***********************************************************************
