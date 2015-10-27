!deck Shell_Info
!***begin prologue     Shell_Info
!***date written       140601   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Shell_Info
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            Read in radial and angular information needed to compute
!***                   the radial points and weights.
!***description  
!***             
!***             
!***references
!***routines called    
!***end prologue       Shell_Info
!********************************************************************************
!********************************************************************************
                        MODULE Shell_Info
  USE Data
  USE Grid_Defined_Types
  USE Matrix_Print
  USE Lebedev_Quadrature
  USE Gauss_Angular_Quadrature
  USE Renormalization
  IMPLICIT NONE
!
!********************************************************************************
!********************************************************************************
!
                            INTERFACE Quadrature_Grids
                       MODULE PROCEDURE Lebedev_Grid,               &
                                        Theta_Grid,                 &  
                                        Phi_Grid,                   &  
                                        Radial_Grid  
                            END INTERFACE Quadrature_Grids
!
                            INTERFACE Generate_Grid_Variables
                       MODULE PROCEDURE Generate_Angular_Variables, &  
                                        Generate_Radial_Variables  
                            END INTERFACE Generate_Grid_Variables
!
                            INTERFACE Set_Grid_Variables
                       MODULE PROCEDURE Set_Angular_Variables,      &  
                                        Set_Radial_Variables  
                            END INTERFACE Set_Grid_Variables
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
  subroutine satshl(atom)
  IMPLICIT NONE
  TYPE(ATOMS)                            :: atom
  INTEGER                                :: lenth 
  INTEGER                                :: intkey
  INTEGER                                :: precision
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: i_sign
  INTEGER                                :: j_sign
  INTEGER                                :: rule
  INTEGER                                :: diff
  INTEGER                                :: len
  CHARACTER (LEN=3)                      :: itoc
  CHARACTER (LEN=80)                     :: chrkey
  LOGICAL                                :: logkey
  LOGICAL                                :: dollar
  REAL(idp)                              :: start
  REAL(idp)                              :: end
  REAL(idp)                              :: increment
  REAL(idp)                              :: fpkey
!
!
! For each center or atom lay out a grid.  
! The grid consists of shells.  In each shell there is an angular and radial set of points.  
! The distribution of both the angular and radial grid is arbitrary.
! The angular grid can be either a non-separable Lebedev grid or a product of theta and phi points.
! If the grid is a lebedev grid then the distribution of points cannot be arbitrary.  
! For a product gris intheta and phi, the breakup into nregions is arbitrary.
  call pakstr(cpass,lstrng)
  IF (cpass(1:lstrng) == 'atom') THEN
      Write(iout,1) str 
  ELSE IF ( cpass(1:lstrng) == 'scattering') THEN
      Write(iout,2) str
  ELSE
      Call lnkerr('error in center pass')
  END IF
!                   
!                           Loop over each atomic center
!
  IF ( dollar('$quadrature_'//str(1:lstrng),card,cpass,inp) ) THEN  ! Get a data card for this center
       CALL cardin(card)
       atom%n_shl = intkey(card,'no_of_shells',1,' ')
       ALLOCATE(atom%shl(1:atom%n_shl))
       nonsep=logkey(card,'angular_quadrature=lebedev',.false.,' ')
       IF (nonsep == .true. ) THEN    ! This is a Lebedev gid in angle
           atom%n_coord=2             ! The angular coordinates are not separable.
           WRITE(iout,3)
           CALL iosys('write character "non separable quadrature "'//str//' to grid',0,0,0,'yes')
           DO i = 1, atom%n_coord - 1
              atom%coord_label='lebedev'     
              call pakstr(atom%coord_label,len)
              DO j = 1, atom%n_shl
                 str='shell'//'_'//itoc(j)
                 call pakstr(str,lstrng)
                 IF ( dollar('$'//atom%coord_label(1:len)//'_shell_'//itoc(j),card,cpass,inp) ) THEN  
                      CALL cardin(card)  ! Get a data card for this coordinate and shell
                 ELSE
                      Call lnkerr('error in processing card')
                 END IF
                 write(iout,5) i, j
                 Call Quadrature_Grids(atom%shl(j)%ang%leb_ang)
              END DO
           END DO
       ELSE                          !  This is a separable grid in the angle theta and phi
           WRITE(iout,4)
           atom%n_coord=3
           CALL iosys('write character "non separable quadrature "'//str//' to grid',0,0,0,'no')
           DO i = 1, atom%n_coord - 1  ! loop over two angular coordinates
              atom%coord_label=chrkey(card,'coordinate_label_'//itoc(i),'theta',' ')     
              call pakstr(atom%coord_label,len)
              IF ( atom%coord_label=='theta') THEN     
                   IF ( dollar('$'//atom%coord_label(1:len),card,cpass,inp) ) THEN  
                        CALL cardin(card)  ! Get a data card for this coordinate and shell
                        n_fixed=intkey(card,'no_of_fixed_points',0,' ')
                        fixed(1)=.false.
                        fixed(2)=.false.
                        drop(1)=.false.
                        drop(2)=.false.
                        IF(n_fixed /= 0) THEN
                           fix(1)=logkey(card,'left_fixed_point',.false.,' ')
                           fix(2)=logkey(card,'right_fixed_point',.false.,' ')
                           drop(1)=logkey(card,'drop_left_function',.false.,' ')
                           drop(2)=logkey(card,'drop_right_function',.false.,' ')
                        END IF
                   ELSE
                        Call lnkerr('error in processing card')
                   END IF
                   DO j = 1, atom%n_shl
                                              ! Loop over the shells
                      str='shell'//'_'//itoc(j)
                      call pakstr(str,lstrng)
                      IF ( dollar('$'//atom%coord_label(1:len)//'_shell_'//itoc(j),card,cpass,inp) ) THEN  
                           CALL cardin(card)  ! Get a data card for this coordinate and shell
                      ELSE
                           Call lnkerr('error in processing card')
                      END IF
                      Call Set_Grid_Variables(atom%shl(j)%ang,atom%shl(j)%ang%theta_ang%reg)
                      write(iout,5) i, j
                      Call Quadrature_Grids(atom%shl(j)%ang,atom%shl(j)%ang%theta_ang,   &
                                            atom%shl(j)%ang%theta_ang%reg)
                   END DO
              ELSE IF( atom%coord_label=='phi') THEN     
                   IF ( dollar('$'//atom%coord_label(1:len),card,cpass,inp) ) THEN  
                        CALL cardin(card)  ! Get a data card for this coordinate and shell
                        n_fixed=intkey(card,'no_of_fixed_points',0,' ')
                        fixed(1)=.false.
                        fixed(2)=.false.
                        drop(1)=.false.
                        drop(2)=.false.
                        IF(n_fixed /= 0) THEN
                           fix(1)=logkey(card,'left_fixed_point',.false.,' ')
                           fix(2)=logkey(card,'right_fixed_point',.false.,' ')
                           drop(1)=logkey(card,'drop_left_function',.false.,' ')
                           drop(2)=logkey(card,'drop_right_function',.false.,' ')
                        END IF
                   ELSE
                        Call lnkerr('error in processing card')
                   END IF
                   DO j = 1, atom%n_shl
!                                                Loop over the shells
                      str='shell'//'_'//itoc(j)
                      call pakstr(str,lstrng)
                      IF ( dollar('$'//atom%coord_label(1:len)//'_shell_'//itoc(j),card,cpass,inp) ) THEN  
                           CALL cardin(card)  ! Get a data card for this coordinate and shell
                      ELSE
                           Call lnkerr('error in processing card')
                      END IF
                      Call Set_Grid_Variables(atom%shl(j)%ang, atom%shl(j)%ang%phi_ang%reg)
                      write(iout,5) i, j
                      Call Quadrature_Grids(atom%shl(j)%ang, atom%shl(j)%ang%phi_ang,                 &
                                            atom%shl(j)%ang%phi_ang%reg)
                   END DO
              ELSE
                   Call lnkerr('bad coordinate card')
              END IF

           END DO
       END IF
!                                                Do the radial coordinate     
       atom%coord_label='radial'
       IF ( dollar('$'//atom%coord_label,card,cpass,inp) ) THEN
            CALL cardin(card)
            n_fixed=intkey(card,'no_of_fixed_points',0,' ')
            fixed(1)=.false.
            fixed(2)=.false.
            drop(1)=.false.
            drop(2)=.false.
            IF(n_fixed /= 0) THEN
               fix(1)=logkey(card,'left_fixed_point',.false.,' ')
               fix(2)=logkey(card,'right_fixed_point',.false.,' ')
               drop(1)=logkey(card,'drop_left_function',.false.,' ')
               drop(2)=logkey(card,'drop_right_function',.false.,' ')
            END IF
       ELSE
            Call lnkerr('error in processing card')
       END IF
       DO j = 1, atom%n_shl
!                                 Loop over the shells 
          str='shell'//'_'//itoc(j)
          call pakstr(str,lstrng)
          IF ( dollar('$'//atom%coord_label(1:len)//'_shell_'//itoc(j),card,cpass,inp) ) THEN
               CALL cardin(card)  ! Get a data card for this coordinate and shell                            
          ELSE
               Call lnkerr('error in processing card')
          END IF
          atom%shl(j)%num=j
          Call Set_Grid_Variables(atom%shl(j)%rad, atom%shl(j)%rad%reg )
          Call Quadrature_Grids( atom%shl(j)%rad,atom%shl(j)%rad%reg)
      END DO
  ELSE
       Call lnkerr('error in quadrature input')
  END IF
!
1  Format(/,25x,' quadrature data for atom = ',i2)
2  Format(/,25x,' quadrature data for scattering center')
3  FORMAT(' angular quadrature is non-separable (lebedev) in angles')
4  FORMAT(' angular quadrature is separable in angles')
5  FORMAT(' quadrature data for coordinate = ',i1,' shell = ',i4)
END SUBROUTINE satshl
!********************************************************************************
!********************************************************************************
!deck Set_Angular_Variables
!***begin prologue     Set_Angular_Variables
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Theta grid
!***author             schneider, barry (nsf)
!***source             
!***purpose            Set boundary conditionsd for first and last shell
!***  
!***routines called
!***end prologue       Set_Angular_Variables
  SUBROUTINE Set_Angular_Variables(ang,reg)
  IMPLICIT NONE
  TYPE(ANGULAR)                             :: ang
  TYPE(REGIONAL), DIMENSION(:), ALLOCATABLE :: reg
  INTEGER                                   :: i
  INTEGER                                   :: j
  INTEGER                                   :: nblock
  CHARACTER (LEN=3)                         :: itoc
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
     ang%n_reg=0
     nblock=intkey(card,'number_of_major_blocks',1,' ')       
     ALLOCATE(work(1:maxreg))
     DO i=1,nblock 
        IF ( dollar('$block_'//itoc(i),card, cpass,inp) ) THEN
             skip=logkey(card,'skip',.false.,' ')
             IF(skip) THEN
                write(iout,1) i
             ELSE 
                nply=intkey(card,'default_order',10,' ')
                nsubr=intkey(card,'number_of_subregions',1,' ')
                boundl=fpkey(card,'left_boundary',0.d0,' ')
                boundr=fpkey(card,'right_boundary',0.d0,' ')
                write(iout,2) i, nply, nsubr, boundl, boundr
                step = ( boundr - boundl ) / nsubr
                DO j=1,nsubr
                   ang%n_reg= ang%n_reg + 1
                   work(ang%n_reg)%edge(1) = boundl
                   work(ang%n_reg)%edge(2) = work(ang%n_reg)%edge(1) + step
                   boundl = work(ang%n_reg)%edge(2)
                   work(ang%n_reg)%n_pts = nply 
                END DO
             END IF
        END IF
     END DO     
     ALLOCATE(reg(1:ang%n_reg))
     DO i = 1, ang%n_reg
        reg(i)%n_pts = work(i)%n_pts
        reg(i)%edge(1:2) = work(i)%edge(1:2)
        reg(i)%type_quadrature = 'lobatto'
        reg(i)%n_fixed = 0 
        reg(i)%first = 1
        reg(i)%last  = reg(i)%n_pts
        reg(i)%n_fun = reg(i)%n_pts
     END DO
  END IF
  DEALLOCATE(work)
!
! Lets fix up the first region
!
! Special case of one region
!
  IF (ang%n_reg == 1 ) THEN
      i = 1   
      IF (drop(1) == .true.) THEN
          reg(i)%first = reg(i)%first + 1
          reg(i)%n_fun = reg(i)%n_fun - 1
      END IF
      IF (drop(2) == .true.) THEN
          reg(i)%last = reg(i)%last - 1
          reg(i)%n_fun = reg(i)%n_fun - 1
      END IF
      IF (n_fixed == 0 ) THEN
         reg(i)%type_quadrature = 'gauss'
         reg(i)%n_fixed = 0
      ELSE IF (n_fixed == 1 ) THEN
         IF ( fixed(1) == .true. ) THEN
              reg(i)%n_fixed = 1
              reg(i)%type_quadrature = 'radau'
         ELSE IF( fixed(2) == .true. ) THEN
              reg(i)%n_fixed = 2
              reg(i)%type_quadrature = 'radau'
         END IF
      ELSE
         reg(i)%type_quadrature = 'lobatto'
      END IF
  ELSE
      i = 1
      IF (drop(1) == .true.) THEN
          reg(i)%first = reg(i)%first + 1
          reg(i)%n_fun = reg(i)%n_fun - 1
      END IF
      IF (fixed(1) == .true. ) THEN
         reg(i)%type_quadrature = 'lobatto'
      ELSE
         reg(i)%type_quadrature = 'radau'
              reg(i)%n_fixed = 2
      END IF
      i = ang%n_reg
      IF (drop(2) == .true.) THEN
          reg(i)%last = reg(i)%last - 1
          reg(i)%n_fun = reg(i)%n_fun - 1
      END IF
      IF (fixed(2) == .true. )THEN
          reg(i)%type_quadrature = 'lobatto'
      ELSE
          reg(i)%type_quadrature = 'radau'
          reg(i)%n_fixed = 1
      END IF
  END IF
1   FORMAT(/,1x,'skipping input block            = ',i4)
2   FORMAT(/,1x, 'block = ',i3,                            &
           /,15x,'quadrature order              = ',i4,    &
           /,15x,'number of subregions          = ',i4,    & 
           /,15x,'left hand boundary            = ',e15.8, &
           /,15x,'right hand boundary           = ',e15.8)
  END SUBROUTINE Set_Angular_Variables
!********************************************************************************
!********************************************************************************
!deck Set_Radial_Variables
!***begin prologue     Set_Radial_Variables
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Theta grid
!***author             schneider, barry (nsf)
!***source             
!***purpose            Set boundary conditionsd for first and last shell
!***  
!***routines called
!***end prologue       Set_Radial_Variables
  SUBROUTINE Set_Radial_Variables(rad,reg)
  IMPLICIT NONE
  TYPE(RADIAL)                                 :: rad
  TYPE(REGIONAL), DIMENSION(:), ALLOCATABLE    :: reg    
  REAL(idp)                                    :: fpkey
  INTEGER                                      :: i
  INTEGER                                      :: j
  CHARACTER (LEN=3)                            :: itoc
  INTEGER                                      :: intkey
  LOGICAL                                      :: logkey
  LOGICAL                                      :: dollar
!
! Set the defaults for all regions.  The defaults are that lobatto quadrature is used, and both the left and right
! points of each element are fixed.  No basis functions are discarded. By enabling lobatto the Gauss routine does
! not even look at the variable n_fixed and simply knows that both endpoints are fiexed
! 
!
  automte=logkey(card,'automate',.false.,' ')
  reuse_sector_information=logkey(card,'reuse_sector_information',.false.,' ')
  IF(.not.automte) THEN
     rad%n_reg=intkey(card,'no_regions',1,' ')
     CALL fparr(card,'radial_boundaries',rad%edge,2,' ')
     DO i = 1, rad%n_reg
        reg(i)%n_pts = intkey(card,'region_'//itoc(i)//'_order',5,' ')
        reg(i)%type_quadrature = 'lobatto'
        reg(i)%n_fixed = 0 
        reg(i)%first = 1
        reg(i)%last  = reg(i)%n_pts
        reg(i)%n_fun = reg(i)%n_pts
     END DO
  ELSE
     nblock=intkey(card,'number_of_major_blocks',1,' ')       
     ALLOCATE(work(1:maxreg))
     DO i=1,nblock 
        IF ( dollar('$block_'//itoc(i),card, cpass,inp) ) THEN
             skip=logkey(card,'skip',.false.,' ')
             IF(skip) THEN
                write(iout,1) i
             ELSE 
                nply=intkey(card,'default_order',10,' ')
                nsubr=intkey(card,'number_of_subregions',1,' ')
                boundl=fpkey(card,'left_boundary',0.d0,' ')
                boundr=fpkey(card,'right_boundary',0.d0,' ')
                write(iout,2) i, nply, nsubr, boundl, boundr
                step = ( boundr - boundl ) / nsubr
                DO j=1,nsubr
                   rad%n_reg = rad%n_reg + 1
                   work(rad%n_reg)%edge(1) = boundl
                   work(rad%n_reg)%edge(2) = work(rad%n_reg)%edge(1) + step
                   boundl = work(rad%n_reg)%edge(2)
                   work(rad%n_reg)%n_pts = nply 
                END DO
             END IF
        END IF
     END DO     
     ALLOCATE(reg(1:rad%n_reg))
     DO i = 1, rad%n_reg
        reg(i)%n_pts = work(i)%n_pts
        reg(i)%edge(1:2) = work(i)%edge(1:2)
        reg(i)%type_quadrature = 'lobatto'
        reg(i)%n_fixed = 0 
        reg(i)%first = 1
        reg(i)%last  = reg(i)%n_pts
        reg(i)%n_fun = reg(i)%n_pts
     END DO
  END IF
  DEALLOCATE(work)
!
! Lets fix up the first region
!
! Special case of one region
!
  IF (rad%n_reg == 1 ) THEN
      i = 1   
      IF (drop(1) == .true.) THEN
          reg(i)%first = reg(i)%first + 1
          reg(i)%n_fun = reg(i)%n_fun - 1
      END IF
      IF (drop(2) == .true.) THEN
          reg(i)%last = reg(i)%last - 1
          reg(i)%n_fun = reg(i)%n_fun - 1
      END IF
      IF (n_fixed == 0 ) THEN
         reg(i)%type_quadrature = 'gauss'
         reg(i)%n_fixed = 0
      ELSE IF (n_fixed == 1 ) THEN
         IF ( fixed(1) == .true. ) THEN
              reg(i)%n_fixed = 1
              reg(i)%type_quadrature = 'radau'
         ELSE IF( fixed(2) == .true. ) THEN
              reg(i)%n_fixed = 2
              reg(i)%type_quadrature = 'radau'
         END IF
      ELSE
         reg(i)%type_quadrature = 'lobatto'
      END IF
  ELSE
      i = 1
      IF (drop(1) == .true.) THEN
          reg(i)%first = reg(i)%first + 1
          reg(i)%n_fun = reg(i)%n_fun - 1
      END IF
      IF (fixed(1) == .true. ) THEN
         reg(i)%type_quadrature = 'lobatto'
      ELSE
         reg(i)%type_quadrature = 'radau'
              reg(i)%n_fixed = 2
      END IF
      i = rad%n_reg
      IF (drop(2) == .true.) THEN
          reg(i)%last = reg(i)%last - 1
          reg(i)%n_fun = reg(i)%n_fun - 1
      END IF
      IF (fixed(2) == .true. )THEN
          reg(i)%type_quadrature = 'lobatto'
      ELSE
          reg(i)%type_quadrature = 'radau'
          reg(i)%n_fixed = 1
      END IF
  END IF
1   FORMAT(/,1x,'skipping input block            = ',i4)
2   FORMAT(/,1x, 'block = ',i3,                            &
           /,15x,'quadrature order              = ',i4,    &
           /,15x,'number of subregions          = ',i4,    & 
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
  INTEGER                                   :: i
  INTEGER                                   :: intkey
  CHARACTER (LEN=30)                        :: chrkey
  CHARACTER (LEN=3)                         :: itoc
  LOGICAL                                   :: avail
!
  i = intkey(card,'lebedev_rule_number',8,' ')
  order=rule_order_table(i)
  avail=.false.
  IF ( rule_logic_table(i) == 1 ) THEN
       avail = .true.
       precision = precision_table (i)
       write (iout,1) rule_order_table(i), precision
       order=rule_order_table(i)
       ang%nleb=precision_table(i)
       ang%nang=order
  ELSE
       Call lnkerr('error in lebedev rule')
  END IF
  WRITE(iout,1) ang%nleb, ang%nang
  CALL iosys ('write integer "lebedev quadrature order '//str//'" to grid',1,ang%nleb,0,' ')
  CALL iosys ('write integer "no lebedev points '//str//'" to grid',1,ang%nang,0,' ')
  atom%maxang=MAX(atom%maxang,ang%nang)
!
  Call  Generate_Lebedev_Points_Weights(ang, leb_rule)
!
1 FORMAT(' Generating Lebedev Quadrature. Lebedev precision = ',i5,2x,   &
         ' Number of lebedev points = ',i5)
  END SUBROUTINE Lebedev_Grid
!********************************************************************************
!********************************************************************************
!deck Theta_Grid.f
!***begin prologue     Theta_Grid
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Theta grid
!***author             schneider, barry (nsf)
!***source             
!***purpose            read in angular information for theta based grid
!***  
!***routines called
!***end prologue       Theta_Grid
  SUBROUTINE Theta_Grid(ang,theta_ang,reg)
  IMPLICIT NONE
  TYPE(ANGULAR)                             :: ang
  TYPE(THETA)                               :: theta_ang
  TYPE(REGIONAL), DIMENSION(:), ALLOCATABLE :: reg
  TYPE(REAL_VECTOR)                         :: type_real_vector
  TYPE(REAL_MATRIX)                         :: type_real_matrix
  INTEGER                                   :: intkey
  INTEGER                                   :: i
  CHARACTER (LEN=30)                        :: chrkey
  CHARACTER (LEN=3)                         :: itoc
!
! We define a single theta grid independent of the m quantum number.  The grid needs to capture the
! polynomic behavior of all P(l,m) on the interval {-1,1} with a factor of ( 1-x*x)^{1/2} removed.  With the factor 
! removed, these functions are all polynomials in {x}.  How well the exact behavior is captured depends on the 
! size of the quadrature.
!
  ang%start = zero
  ang%end = pi
  ang%increment = ang%end / ang%n_reg
  Call Generate_Grid_Variables(ang,reg)
  DO i= 1, ang%n_reg 
     CALL iosys('write character "no theta quadrature points '//str//' region '//itoc(i)//'" to grid',  &
                 0,0,0,reg(i)%n_pts)
     CALL iosys('write character "theta quadrature type '//str//' region '//itoc(i)//'" to grid',0,0,0, &
                 reg(i)%type_quadrature)
     Write(iout,1) i, reg(i)%n_pts, reg(i)%type_quadrature
  END DO
!
! Calculate factors needed for matrix elements
!
  DO i = 1, ang%n_reg
     ALLOCATE( reg(i)%q_fac(1:reg(i)%n_pts), reg(i)%inv_q_fac(1:reg(i)%n_pts),                           &
               reg(i)%inv_sqrt_q_fac(1:reg(i)%n_pts))
     reg(i)%q_fac(:) = one - reg(i)%q(:) * reg(i)%q(:)
     reg(i)%inv_q_fac(:) = one / reg(i)%q_fac(:)
     reg(i)%inv_sqrt_q_fac(:) = Sqrt ( reg(i)%inv_q_fac(:) )
  END DO
!
  IF (prnt(3) == .true.) THEN
      DO i = 1, ang%n_reg
         call Print_Matrix(type_real_vector,reg(i)%q_fac,title='q_factor Region-'//itoc(i))
         call Print_Matrix(type_real_vector,reg(i)%inv_q_fac,title='inverse_q_factor Region-'//itoc(i))
         call Print_Matrix(type_real_vector,reg(i)%inv_sqrt_q_fac,title='inverse_sqrt_q_factor Region-'//itoc(i))
     END DO
!
  END IF  
!
! Get normalization factors and renormalize
!
  Call Normalization(reg,ang%n_reg)
  Call Normalized_Functions(reg,ang%n_reg)
!
  IF (prnt(2) == .true.) THEN
      DO i = 1, ang%n_reg
         call Print_Matrix(type_real_matrix,reg(i)%bf(:,reg(i)%first:reg(i)%last),                  &
                           reg(i)%n_pts,reg(i)%n_fun,                                               &
                           title='Final Normalized Basis Polynomials Region-'//itoc(i))
         call Print_Matrix(type_real_matrix,reg(i)%dbf(:,reg(i)%first:reg(i)%last),                 &
                           reg(i)%n_pts,reg(i)%n_fun,                                               &
                           title='First Derivative of Final Normalized Basis Polynomials Region-'//itoc(i))
         call Print_Matrix(type_real_matrix,reg(i)%ddbf(:,reg(i)%first:reg(i)%last),                &
                           reg(i)%n_pts,reg(i)%n_fun,                                               &
                           title='Second Derivative of Final Normalized Basis Polynomials Region-'//itoc(i))
     END DO
  END IF  
1 FORMAT(' theta region = ',i3,2x,' number of points = ',i5,2x,' type quadrature = ',a8)
  END SUBROUTINE Theta_Grid
!********************************************************************************
!********************************************************************************
!deck Phi_Grid.f
!***begin prologue     Phi_Grid
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Phi grid
!***author             schneider, barry (nsf)
!***source             
!***purpose            read in angular information for theta based grid
!***  
!***routines called
!***end prologue       Phi_Grid
  SUBROUTINE Phi_Grid(ang,phi_ang,reg)
  IMPLICIT NONE
  TYPE(ANGULAR)                             :: ang
  TYPE(PHI)                                 :: phi_ang
  TYPE(REGIONAL), DIMENSION(:), ALLOCATABLE :: reg
  TYPE(REAL_VECTOR)                         :: type_real_vector
  TYPE(REAL_MATRIX)                         :: type_real_matrix
  INTEGER                                   :: intkey
  INTEGER                                   :: i
  INTEGER                                   :: ndef
  CHARACTER (LEN=30)                        :: chrkey
  CHARACTER (LEN=3)                         :: itoc
!
  ang%start = 0.d0  
  ang%end = two_pi
  ang%increment = ang%end / ang%n_reg
  Call Generate_Grid_Variables(ang,reg)
  DO i= 1, ang%n_reg 
     CALL iosys('write character "phi quadrature order '//str//' region '//itoc(i)//'" to grid',0,0,0,          &
                 reg(i)%n_pts)
     CALL iosys('write character "phi quadrature type '//str//' region '//itoc(i)//'" to grid',0,0,0,           &
                 reg(i)%type_quadrature)
     Write(iout,1) i, reg(i)%n_pts, reg(i)%type_quadrature
  END DO
!
  DO i = 1, ang%n_reg
     ALLOCATE( reg(i)%q_fac(1:reg(i)%n_pts), reg(i)%inv_q_fac(1:reg(i)%n_pts),                           &
               reg(i)%inv_sqrt_q_fac(1:reg(i)%n_pts))
     reg(i)%q_fac(:) = one
     reg(i)%inv_q_fac(:) = one
     reg(i)%inv_sqrt_q_fac(:) = one
  END DO
  IF (prnt(3) == .true.) THEN
      DO i = 1, ang%n_reg
         call Print_Matrix(type_real_vector,reg(i)%q_fac,title='q_factor Region-'//itoc(i))
         call Print_Matrix(type_real_vector,reg(i)%inv_q_fac,title='inverse_q_factor Region-'//itoc(i))
         call Print_Matrix(type_real_vector,reg(i)%inv_sqrt_q_fac,title='inverse_sqrt_q_factor Region-'//itoc(i))
     END DO
  END IF  
!
! Get normalization factors and renormalize
!
  Call Normalization(reg,ang%n_reg)
  Call Normalized_Functions(reg,ang%n_reg)
!
  IF (prnt(2) == .true.) THEN
      DO i = 1, ang%n_reg
         call Print_Matrix(type_real_matrix,reg(i)%bf(:,reg(i)%first:reg(i)%last),                  &
                           reg(i)%n_pts,reg(i)%n_fun,                                               &
                           title='Final Normalized Basis Polynomials Region-'//itoc(i))
         call Print_Matrix(type_real_matrix,reg(i)%dbf(:,reg(i)%first:reg(i)%last),                 &
                           reg(i)%n_pts,reg(i)%n_fun,                                               &
                           title='First Derivative of Final Normalized Basis Polynomials Region-'//itoc(i))
         call Print_Matrix(type_real_matrix,reg(i)%ddbf(:,reg(i)%first:reg(i)%last),                &
                           reg(i)%n_pts,reg(i)%n_fun,                                               &
                           title='Second Derivative of Final Normalized Basis Polynomials Region-'//itoc(i))
     END DO
  END IF  
1 FORMAT(' phi region = ',i3,2x,' number of points = ',i5,2x,' type quadrature = ',a8)
  END SUBROUTINE Phi_Grid
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
  SUBROUTINE Radial_Grid(rad,reg)
  IMPLICIT NONE
  TYPE(RADIAL)                                 :: rad
  TYPE(REGIONAL), DIMENSION(:), ALLOCATABLE    :: reg
  TYPE(REAL_VECTOR)                            :: type_real_vector
  TYPE(REAL_MATRIX)                            :: type_real_matrix
  INTEGER                                      :: intkey
  INTEGER                                      :: i
  LOGICAL                                      :: dollar
  CHARACTER (LEN=30)                           :: chrkey
  CHARACTER (LEN=3)                            :: itoc
!
  DO i = 1, rad%n_reg
     str='shell'//'_'//itoc(i)
     call pakstr(str,lstrng)
     IF ( dollar('$'//atom%coord_label//'_region_'//itoc(i),card,cpass,inp) ) THEN  
          CALL cardin(card)  ! Get a data card for this coordinate and shell
     ELSE
          Call lnkerr('error in processing card')
     END IF
     Call Generate_Grid_Variables(rad,reg(i)) 
     CALL iosys('write character "no radial quadrature points '//str//' to grid',0,0,0,   &
                 reg(i)%n_pts)
     CALL iosys('write character "theta quadrature type '//str//'" to grid',0,0,0,        &
                 reg(i)%type_quadrature)
     Write(iout,1) reg(i)%n_pts, reg(i)%type_quadrature
  END DO
!
! Coordinate factors needed for matrix elements
!
  DO i = 1, rad%n_reg
     ALLOCATE( reg(i)%q_fac(1:reg(i)%n_pts), reg(i)%inv_q_fac(1:reg(i)%n_pts),            &
               reg(i)%inv_sqrt_q_fac(1:reg(i)%n_pts))
     reg(i)%q_fac(:) = reg(i)%q(:) * reg(i)%q(:)
     reg(i)%inv_q_fac(:) = one / reg(i)%q_fac(:)
     reg(i)%inv_sqrt_q_fac(:) = one / reg(i)%q(:)
  END DO
!
  IF (prnt(3) == .true.) THEN
      DO i = 1, rad%n_reg
         call Print_Matrix(type_real_vector,reg(i)%q_fac,title='q_factor'//str)
         call Print_Matrix(type_real_vector,reg(i)%inv_q_fac,title='inverse_q_factor'//str)
         call Print_Matrix(type_real_vector,reg(i)%inv_sqrt_q_fac,title='inverse_sqrt_q_factor'//str)
     END DO
  END IF  
!
! Get normalization factors and renormalize
!
  Call Normalization(reg,rad%n_reg)
  Call Normalized_Functions(reg,rad%n_reg)
!
1 FORMAT(' radial region = ',i3,2x,' number of points = ',i5,2x,' type quadrature = ',a8)
  END SUBROUTINE Radial_Grid
!********************************************************************************
!********************************************************************************
!deck Generate_Angular_Variables
!***begin prologue     Generate_Angular_Variables
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Angular
!***author             schneider, barry (nsf)
!***source             
!***purpose            Generate_Angular_Variables
!***routines called
!***end prologue       Generate_Angular_Variables
  SUBROUTINE Generate_Angular_Variables(ang,reg)
  IMPLICIT NONE
  TYPE(ANGULAR)                             :: ang
  TYPE(REGIONAL), DIMENSION(:), ALLOCATABLE :: reg
  TYPE(REAL_VECTOR)                         :: type_real_vector
  TYPE(REAL_MATRIX)                         :: type_real_matrix
  INTEGER                                   :: i
  INTEGER                                   :: intkey
  INTEGER                                   :: n_reg
  CHARACTER (LEN=3)                         :: itoc
!
  DO i= 1, ang%n_reg 
     reg(i)%edge(1) = ang%start
     reg(i)%edge(2) = ang%end
     Call fparr (card,'edges_region_'//itoc(i),reg(i)%edge,2,' ')
     ang%start = ang%end
     ang%end = ang%end + ang%increment
  END DO
!
  DO i= 1, ang%n_reg 
     ALLOCATE( reg(i)%q(1:reg(i)%n_pts), reg(i)%wt(1:reg(i)%n_pts) )
     Call gauss(reg(i)%q, reg(i)%wt, reg(i)%edge,type_quadrature=reg(i)%type_quadrature,            &
                fixed_point=reg(i)%n_fixed, n=reg(i)%n_pts)
  END DO
  IF (prnt(1) == .true.) THEN
      DO i = 1, ang%n_reg
         call Print_Matrix(type_real_vector,reg(i)%q,title='Points Region-'//itoc(i))
         call Print_Matrix(type_real_vector,reg(i)%wt,title='Weights Region-'//itoc(i))
      END DO
  END IF
  DO i= 1, ang%n_reg 
     ALLOCATE( reg(i)%p(1:reg(i)%n_pts,1:reg(i)%n_pts),reg(i)%dp(1:reg(i)%n_pts,1:reg(i)%n_pts),    &
               reg(i)%ddp(1:reg(i)%n_pts,1:reg(i)%n_pts),reg(i)%normalization(1:reg(i)%n_pts) )
     CALL cpoly(reg(i)%p,reg(i)%dp,reg(i)%ddp,reg(i)%q,reg(i)%n_pts)
  END DO
  IF (prnt(2) == .true.) THEN
      DO i = 1, n_reg
         call Print_Matrix(type_real_matrix,reg(i)%p,reg(i)%n_pts,reg(i)%n_pts,                     &
                           title='Unnormalized Polynomials Region-'//itoc(i))
         call Print_Matrix(type_real_matrix,reg(i)%dp,reg(i)%n_pts,reg(i)%n_pts,                    &
                           title='First Derivative of Unnormalized Polynomials Region-'//itoc(i))
         call Print_Matrix(type_real_matrix,reg(i)%ddp,reg(i)%n_pts,reg(i)%n_pts,                   &
                           title='Second Derivative of Unnormalized Polynomials Region-'//itoc(i))
     END DO
  END IF                                                                                        
!
  END SUBROUTINE Generate_Angular_Variables
!********************************************************************************
!********************************************************************************
!deck Generate_Radial_Variables
!***begin prologue     Generate_Radial_Variables
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Radial_Radial, link m6200
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            Generate_Radial_Variables
!***routines called
!***end prologue       Generate_Radial_Variables
  SUBROUTINE Generate_Radial_Variables(rad,reg)
  IMPLICIT NONE
  TYPE(RADIAL)                              :: rad
  TYPE(REGIONAL)                            :: reg
  TYPE(REAL_VECTOR)                         :: type_real_vector
  TYPE(REAL_MATRIX)                         :: type_real_matrix
  REAL(idp)                                 :: start
  REAL(idp)                                 :: end
  INTEGER                                   :: i
  INTEGER                                   :: intkey
  CHARACTER (LEN=3)                         :: itoc

!
  Call gauss(reg%q, reg%wt, reg%edge, type_quadrature=reg%type_quadrature,              &
             fixed_point=reg%n_fixed, n=reg%n_pts)
  IF (prnt(1) == .true.) THEN
      call Print_Matrix(type_real_vector,reg%q,title='Points Region-'//itoc(i))
      call Print_Matrix(type_real_vector,reg%wt,title='Weights Region-'//itoc(i))
  END IF
  ALLOCATE( reg%p(1:reg%n_pts,1:reg%n_pts), reg%dp(1:reg%n_pts,1:reg%n_pts),            &
            reg%ddp(1:reg%n_pts,1:reg%n_pts),reg%normalization(1:reg%n_pts) )
!
  CALL cpoly(reg%p,reg%dp,reg%ddp,reg%q,reg%n_pts)
!                  
  IF (prnt(2) == .true.) THEN
      call Print_Matrix(type_real_matrix,reg%p,reg%n_pts,reg%n_pts,                     &
                        title='Unnormalized Polynomials Region-'//itoc(i))
      call Print_Matrix(type_real_matrix,reg%dp,reg%n_pts,reg%n_pts,                    &
                        title='First Derivative of Unnormalized Polynomials Region-'//itoc(i))
      call Print_Matrix(type_real_matrix,reg%ddp,reg%n_pts,reg%n_pts,                   &
                        title='Second Derivative of Unnormalized Polynomials Region-'//itoc(i))
  END IF                                                                                                             
!
  END SUBROUTINE Generate_Radial_Variables
!***********************************************************************
!***********************************************************************
END MODULE Shell_Info
