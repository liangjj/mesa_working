!***********************************************************************
! Read_Grid_Data
!**begin prologue     Read_Grid_Data
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           
!***
!**author             schneider, b. i.(nsf)
!**source             
!**purpose            Read in Grid Data
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      Read_Grid_Data
!***********************************************************************
!***********************************************************************
                           MODULE Read_Grid_Data
                           USE Data
!                           USE FEDVR_Shared
                           USE Grid_Defined_Types
!                           USE Matrix_Print
                    IMPLICIT NONE
!
!***********************************************************************
!***********************************************************************                          
!                            Explicit Interfaces

                            INTERFACE Read_Coordinate_Information
                       MODULE PROCEDURE Read_Cartesian,                    &
                                        Read_Spherical,                    &
                                        Read_Cylindrical,                  &
                                        Read_Spheroidal
                            END INTERFACE Read_Coordinate_Information
!***********************************************************************
!***********************************************************************
                           Contains
!***********************************************************************
!***********************************************************************
!deck Set_General_Parameters
!***begin prologue     Set_General_Parameters
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input routines for grid and dvr basis sets
!***description        Note that the various Xkey functions enable keyword
!                      input, with default values if the keyword does not appear 
!                      in the input line.  There are Xkey routine for real, integer,
!                      logical and character type variables.
!                      The general process is that the card variable is read and
!                      the string decoded by the Xkey routine which set the variables.
!
!                      The prnkey parameter is;
!                      'sector_points', 'sector_factors', 'sector_polynomials', 
!                      'sector_matrices' ,'global_points', 'global_polynomials', 
!                      'potential','global_matrices', 'hamiltonian', 'eigenvalues',
!                      'eigenvectors', 'all'
!***references
!***routines called    
!***end prologue       Set_General_Parameters
  SUBROUTINE Set_General_Parameters  
  IMPLICIT NONE
  REAL(idp)                   :: fpkey
  LOGICAL                     :: dollar
  LOGICAL                     :: logkey
  INTEGER                     :: intkey
  CHARACTER(LEN=80)           :: chrkey
  IF ( dollar('$general_keywords',card,cpass,inp) )THEN
       niter=intkey(ops,'grid=number-of-voronoi-iterations',3,' ')
       alpha=fpkey(ops,'grid=exponential-scattering-grid-cutoff', 1.d0,' ')
       prnt=logkey(ops,'grid=print',.false.,' ')
       no_voroni=logkey(ops,'grid=voronoi=off',.false.,' ')
       cutoff=fpkey(ops,'grid=scattering-grid-cutoff',10.d0,' ')
       no_scat=logkey(ops,'grid=no-scattering-center',.false.,' ')  ! Default is scattering center present                  
       one_grid=logkey(ops,'grid=only-one-grid',.false.,' ')        ! Default is multiple grids                           
       IF (one_grid) THEN  ! If only one grid, that grid is the scattering grid centered at (0,.0.,0.)                       
           no_scat=.false.
           no_voroni=.true. ! No Voronoi transformation                                                                      
       END IF
       no_disk=logkey(ops,'grid=no-disk-output',.false.,' ')
       yukawa_on=logkey(ops,'grid=yukawa=on',.false.,' ')
       yn='yes'
       IF (no_scat==.true.) THEN
           yn='true'
       END IF
  ELSE
       call lnkerr('General Keywords Absent:Quit')
  END If
  Write(iout,1) 
  CALL iosys ('read character "grid filename" from rwf',-1,0,0,grid_filename) ! Open the grid file                            
  CALL iosys ('open grid as new',0,0,0,grid_filename)
  CALL iosys('write character "scattering center" to grid',0,0, 0,yn)  ! Start writing data to the grid file                 
!            read centers and the parameters for the yukawa potential                                                         
  IF ( dollar('$centers',card,cpass,inp) ) THEN
       CALL cardin(card)
  ELSE
       Call lnkerr('error in centers input')
  END IF
  ncent=intkey(card,'no-atomic-centers',1,' ')
  ncplus=ncent
  IF (no_scat==.false.) THEN
      ncplus = ncent + 1
  END IF
  CALL iosys ('write integer "number of atomic centers" to grid',1,ncent,0,' ')
  CALL iosys ('create real "atomic center positions" on grid',3*ncent,0,0,' ')
  CALL iosys ('create real "yukawa exponents" on grid',ncent,0,0,' ')
  CALL iosys ('create real "nuclear charges" on grid',ncent,0,0,' ')
  ALLOCATE(atom(1:ncplus))
  DO i=1,ncent
     ALLOCATE(atom(i)%a(1:3))
     atom(i)%eta=fpkey(card,'exponent-center-'//itoc(i),0.d+00,' ')
     atom(i)%znuc=fpkey(card,'charge-center-'//itoc(i),1.d+00,' ')
     CALL fparr(card,'position-center-'//itoc(i),atom(i)%a, 3,' ')
     WRITE(iout,2) atom(i)%eta, atom(i)%a(1:3), atom(i)%znuc
     CALL iosys ('write real "atomic atomer positions" to grid without rewinding',3,atom(i)%a,0,' ')
     CALL iosys ('write real "yukawa exponents" to grid without rewinding',1,atom(i)%eta,0,' ')
     CALL iosys ('write real "nuclear charges" to grid without rewinding',1,atom(i)%znuc,0,' ')
  END DO
!          scattering center coordinates are read in here. the default is (0.,0.,0.)                                         
  IF (no_scat==.false.) THEN
      CALL fparr(card,'scattering-center-position',atom(ncplus)%a,3,' ')
      WRITE(iout,3) atom(ncplus)%a(1:3)
      CALL iosys ('write real "scattering center position" to grid',3,atom(ncplus)%a,0,' ')
  END IF
1 FORMAT (//,15X,'***** grid generation program *****',//)
2 FORMAT(' yukawa exponent = ',f10.5,' center = (',f10.5,',',f10.5,  &
         ',',f10.5')'/,' charge = ',f10.5)
3 FORMAT(' scattering center ',3F10.5)
4 FORMAT('total number of grid points',1X,i5)
5 FORMAT(/,'maximum radial points in a shell ',i4,  &
        1X,'maximum theta points ',i4,/, 'maximum phi points ',i4)
  END SUBROUTINE Set_General_Keywords
!***********************************************************************
!***********************************************************************e
!deck Read_Data
!***begin prologue     Read_Data
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       Read_Data
!  SUBROUTINE Read_Data(grid)
  SUBROUTINE Read_Data
  IMPLICIT NONE
!
  IF(coordinate_system == 'cartesian') THEN
     Call Read_Coordinate_Information(grid%xyz)
  ELSE IF(coordinate_system == 'spherical') THEN
     Call Read_Coordinate_Information(grid%r_theta)
  ELSE IF(coordinate_system == 'cylindrical') THEN
     Call Read_Coordinate_Information(grid%rho_z)
  ELSE IF(coordinate_system == 'spheroidal') THEN
     Call Read_Coordinate_Information(grid%xi_eta)
  END IF  
  END SUBROUTINE Read_Data
!***********************************************************************
!***********************************************************************
!deck Read_Grid_Parameters
!***begin prologue     Read_Grid_Parameters
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       Read_Grid_Parameters
  SUBROUTINE Read_Grid_Parameters
  IMPLICIT NONE
  CHARACTER(LEN=3)                       :: itoc
  INTEGER                                :: intkey
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: nply
  INTEGER                                :: nsubr
  INTEGER                                :: nblock
  LOGICAL                                :: logkey
  LOGICAL                                :: automte
  LOGICAL                                :: dollar
  REAL(idp)                              :: fpkey
  REAL(idp)                              :: step
  REAL(idp)                              :: boundl
  REAL(idp)                              :: boundr
  automte=logkey(card,'automate',.false.,' ')
  reuse_sector_information=logkey(card,'reuse_sector_information',.false.,' ')
  grid%n_fix=intkey(card,'number_of_fixed_points',0,' ')
  grid%fix(1)=.false.
  grid%fix(2)=.false.
  grid%drop(1)=.false.
  grid%drop(2)=.false.
  IF(grid%n_fix /= 0) THEN
     grid%fix(1)=logkey(card,'left_fixed_point',.false.,' ')
     grid%fix(2)=logkey(card,'right_fixed_point',.false.,' ')
     grid%drop(1)=logkey(card,'drop_left_function',.false.,' ')
     grid%drop(2)=logkey(card,'drop_right_function',.false.,' ')
  END IF
  IF(.not.automte) THEN
     nreg=intkey(card,'number_of_regions',1,' ')
     ALLOCATE(grid%reg(1:nreg))
     DO i = 1, nreg
        ALLOCATE( grid%reg(i)%edge(1:2) )
        CALL fparr(card,'region_'//itoc(i)//'_boundaries',grid%reg(i)%edge,2,' ')
        grid%reg(i)%n_pts = intkey(card,'region_'//itoc(i)//'_order',5,' ')
     END DO
  ELSE
     write(iout,*)
     write(iout,*)
     write(iout,*) '                   Automated Selection of', &
                   ' Steps   '
     nreg=0
     nblock=intkey(card,'number_of_major_blocks',1,' ')       
     ALLOCATE(grid%temp(1:maxreg))
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
                   nreg = nreg + 1
                   ALLOCATE(grid%temp(nreg)%edge(1:2))
                   grid%temp(nreg)%edge(1) = boundl
                   grid%temp(nreg)%edge(2) = grid%temp(nreg)%edge(1) + step
                   boundl = grid%temp(nreg)%edge(2)
                   grid%temp(nreg)%n_pts = nply 
                END DO
             END IF
        END IF
     END DO     
     ALLOCATE(grid%reg(1:nreg))
     DO i = 1, nreg
        ALLOCATE(grid%reg(i)%edge(1:2))
        grid%reg(i)%n_pts = grid%temp(i)%n_pts
        grid%reg(i)%edge(1:2) = grid%temp(i)%edge(1:2)
        grid%reg(i)%type_quadrature = 'lobatto'
        grid%reg(i)%n_fixed = int_two
     END DO
     IF ( nreg == int_one ) THEN
          IF ( grid%n_fix == 0 ) THEN
               grid%reg(1)%type_quadrature = 'gauss'
               grid%reg(i)%n_fixed = int_zero
          ELSE IF ( grid%n_fix == 1 ) THEN
               grid%reg(1)%type_quadrature = 'radau'
               grid%reg(i)%n_fixed = int_one
               IF ( grid%fix(1) == .true. ) THEN
                    grid%reg(1)%fixed_point = int_one
               ENDIF
               IF ( grid%fix(2) == .true. ) THEN
                    grid%reg(1)%fixed_point = int_two
               END IF
          ELSE IF  ( grid%n_fix == 2 ) THEN
               grid%reg(1)%type_quadrature = 'lobatto'
          END IF
     ELSE
          grid%reg(1)%type_quadrature    = 'radau'
          grid%reg(1)%fixed_point = int_two
          grid%reg(nreg)%type_quadrature = 'radau'
          grid%reg(nreg)%fixed_point = int_one
          IF ( grid%fix(1) == .true. ) THEN
               grid%reg(1)%type_quadrature = 'lobatto'
               grid%reg(1)%fixed_point = int_one
          ENDIF
          IF ( grid%fix(2) == .true. ) THEN
               grid%reg(nreg)%type_quadrature = 'lobatto'
               grid%reg(nreg)%fixed_point = int_two
          END IF
     END IF          
!
     DEALLOCATE(grid%temp)
     Write(iout,3)
     DO i = 1, nreg
        write(iout,4) i, grid%reg(i)%n_pts, grid%reg(i)%edge(1:2),  &
                      grid%reg(i)%n_fixed,  grid%reg(i)%type_quadrature
     END DO
  END IF
1   FORMAT(/,1x,'skipping input block            = ',i4)
2   FORMAT(/,1x, 'block = ',i3,                            &
           /,15x,'quadrature order              = ',i4,    &
           /,15x,'number of subregions          = ',i4,    & 
           /,15x,'left hand boundary            = ',e15.8, &
           /,15x,'right hand boundary           = ',e15.8)
3   FORMAT(/,1x,'region',5x,'number of points',5x,'left edge',5x,'right edge', &
             7x,'fixed points',5x,'type quadrature' )
4 Format(i4,12x,i4,10x,f12.8,3x,f12.8,10x,i2,17x,a8)
  END SUBROUTINE Read_Grid_Parameters
!***********************************************************************
!***********************************************************************
!deck Read_Cartesian
!***begin prologue     Read_Cartesian
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            
!***description        
!***                   
!***references
!***routines called    
!***end prologue       Read_Cartesian
  SUBROUTINE Read_Cartesian(xyz)
  IMPLICIT NONE
  TYPE(cartesian)             :: xyz
  LOGICAL                     :: dollar
  LOGICAL                     :: logkey
  INTEGER                     :: intkey
  CHARACTER(LEN=80)           :: chrkey
  IF ( dollar('$cartesian',card,cpass,inp) ) THEN
       reuse=logkey(card,'reuse_space_data',.false.,' ')
       xyz%axis = chrkey(card,'axis','x',' ')
!
       Call Read_Grid_Parameters
  ELSE
       call lnkerr('no coordinate keyword found for label = cartesian')
  END If
  write(iout,1) grid%n_fix, grid%fix, grid%drop
!
1 Format(/,15x,'Cartesian Coordinates',15x,                              &
               'number of fixed points    = ',i1,/15x,                   &
               'left point fixed          = ',l1,/,15x,                  &
               'right point fixed         = ',l1,/,15x,                  &
               'left point dropped        = ',l1,/15x,                   &
               'right point dropped       = ',l1)
  END SUBROUTINE Read_Cartesian
!***********************************************************************
!***********************************************************************
!deck Read_Spherical
!***begin prologue     Read_Spherical
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            
!***description        
!***                   
!***references
!***routines called    
!***end prologue       Read_Spherical
  SUBROUTINE Read_Spherical(r_theta)
  IMPLICIT NONE
  TYPE(spherical)             :: r_theta
  LOGICAL                     :: dollar
  LOGICAL                     :: logkey
  INTEGER                     :: intkey
  CHARACTER(LEN=80)           :: chrkey
  IF ( dollar('$spherical',card,cpass,inp) ) THEN
       r_theta%axis = chrkey(card,'axis','r',' ')
       len=lenth(r_theta%axis)
       r_theta%l_max=intkey(card,'maximum_l_value',0,' ')
       r_theta%m_max=intkey(card,'maximum_m_value',0,' ')
       reuse=logkey(card,'reuse_space_data',.false.,' ')
       drctv=chrkey(card,'type_calculation','all_integrals',' ')
!
       Call Read_Grid_Parameters
!
       IF (r_theta%axis(1:len) == 'r') THEN
           write(iout,1) r_theta%axis(1:len), grid%n_fix, grid%fix, grid%drop, r_theta%l_max
       ELSE IF (r_theta%axis(1:len) == 'theta') THEN
           write(iout,2) r_theta%axis(1:len), grid%n_fix, grid%fix, grid%drop, r_theta%m_max
       END IF
!
  ELSE
       call lnkerr('no coordinate keyword found for label = spherical')
  END IF
1 Format(/,15x,'spherical coordinates:Coordinate = ',a8,                        &
               'number of fixed points           = ',i1,/15x,                   &
               'left point fixed                 = ',l1,/,15x,                  &
               'right point fixed                = ',l1,/,15x,                  &
               'left point dropped               = ',l1,/15x,                   &
               'right point dropped              = ',l1,/15x,                   &
               'maximum l value                  = ',i4)
2 Format(/,15x,'spherical coordinates:Coordinate = ',a8,                        &
               'number of fixed points           = ',i1,/15x,                   &
               'left point fixed                 = ',l1,/,15x,                  &
               'right point fixed                = ',l1,/,15x,                  &
               'left point dropped               = ',l1,/15x,                   &
               'right point dropped              = ',l1,/15x,                   &
               'maximum m value                  = ',i4) 
  END SUBROUTINE Read_Spherical
!***********************************************************************
!***********************************************************************
!deck Read_Cylindrical
!***begin prologue     Read_Cylindrical
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            
!***description        
!***                   
!***references
!***routines called    
!***end prologue       Read_Cylindrical
  SUBROUTINE Read_Cylindrical(rho_z)
  IMPLICIT NONE
  TYPE(cylindrical)           :: rho_z
  LOGICAL                     :: dollar
  LOGICAL                     :: logkey
  INTEGER                     :: intkey
  CHARACTER(LEN=80)           :: chrkey
!
  IF ( dollar('$cylindrical',card,cpass,inp) ) THEN
       rho_z%axis = chrkey(card,'axis','rho',' ')
       len=lenth(rho_z%axis)
       reuse=logkey(card,'reuse_space_data',.false.,' ')
       rho_z%m_max=intkey(card,'maximum_m_value',0,' ')
       Call Read_Grid_Parameters
!
!      Determine nphy
!
       write(iout,1) grid%n_fix, grid%fix, grid%drop, rho_z%m_max
!
  ELSE
       call lnkerr('no coordinate keyword found for label = cylindrical')
  END If
1 Format(/,15x,'cylindrical coordinates:Coordinate = ',a8,                        &
               'number of fixed points             = ',i1,/15x,                   &
               'left point fixed                   = ',l1,/,15x,                  &
               'right point fixed                  = ',l1,/,15x,                  &
               'left point dropped                 = ',l1,/15x,                   &
               'right point dropped                = ',l1,/15x,                   &
               'maximum m value                    = ',i4)  
  END SUBROUTINE Read_Cylindrical
!*************************************************************************************
!*************************************************************************************
!deck Read_Spheroidal
!***begin prologue     Read_Spheroidal
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for spheroidal dvr basis sets
!***description        sets up the unperturbed dv hamiltonian for a two center
!***                   problem where charge Z_a is at -R/2 and Z_b at +R/2
!***references
!***routines called    
!***end prologue       Read_Spheroidal
  SUBROUTINE Read_Spheroidal(xi_eta)
  IMPLICIT NONE
  TYPE(spheroidal)                 :: xi_eta
  LOGICAL                          :: dollar
  REAL(idp)                        :: fpkey
  INTEGER                          :: intkey
  CHARACTER(LEN=80)                :: chrkey
!
  IF ( dollar('$spheroidal')) THEN
       xi_eta%axis = chrkey(card,'axis','xi',' ')
       xi_eta%m_max=intkey(card,'maximum_m_value',0,' ')
       xi_eta%Z_a = fpkey(card,'nuclear_charge_left_nucleus',-1.0,' ')
       xi_eta%Z_b = fpkey(card,'nuclear_charge_right_nucleus',1.0,' ')
       xi_eta%R_ab = fpkey(card,'internuclear_separation',2.0,' ')
!
       Call Read_Grid_Parameters
!
!
  ELSE
       call lnkerr('no coordinate keyword found for label = spheroidal')
  END IF
  write(iout,1) grid%n_fix, grid%fix, grid%drop, xi_eta%m_max,           &
                xi_eta%Z_a, xi_eta%Z_b, xi_eta%R_ab
1 Format(/,15x,'spheroidal coordinates:Coordinate = ',a8,                        &
               'number of fixed points            = ',i1,/15x,                   &
               'left point fixed                  = ',l1,/,15x,                  &
               'right point fixed                 = ',l1,/,15x,                  &
               'left point dropped                = ',l1,/15x,                   &
               'right point dropped               = ',l1,/15x,                   &
               'maximum m value                   = ',i4,/15x,                   &
               'charge on nucleus A               = ',f10.5,/15x,                & 
               'charge on nucleus b               = ',f10.5,/15x,                &
               'internuclear distance             = ',e15.8)
  END SUBROUTINE Read_Spheroidal

!***********************************************************************
!***********************************************************************
           END MODULE Read_Grid_Data
!***********************************************************************
!***********************************************************************
