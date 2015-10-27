!deck   m6203
!**begin prologue     m6203
!**date written       
!**revision date      yymmdd   (yymmdd)
!**keywords           
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**description        
!**                   
!**                   
!**references
!**routines called    
!                     
!
  PROGRAM m6203
  USE Lebedev_Quadrature
  USE Gauss_Angular_Quadrature
  USE Grid_Defined_Types
  USE Grid_Generation
  IMPLICIT NONE
  CHARACTER (LEN=3)                        :: itoc
  CHARACTER (LEN=80)                       :: chrkey
  LOGICAL                                  :: dollar
  LOGICAL                                  :: logkey
  LOGICAL                                  :: avail
  LOGICAL                                  :: test_all
  LOGICAL                                  :: type
  REAL(idp)                                :: fpkey
  REAL(idp), DIMENSION(4)                  :: range
  INTEGER                                  :: intkey
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: iloc
  INTEGER                                  :: st
  INTEGER                                  :: irange
  INTEGER                                  :: len
  INTEGER, DIMENSION(:), ALLOCATABLE       :: rule_list
  INTEGER, DIMENSION(:), ALLOCATABLE       :: rule_test
  INTEGER                                  :: number_of_rules_to_test
  CHARACTER (LEN = 80)                     :: title
  CHARACTER (LEN = 80)                     :: type_quadrature
  TYPE(LEBEDEV)                            :: leb_ang
  TYPE(THETA)                              :: theta_ang
  TYPE(PHI)                                :: phi_ang
  TYPE(Integer_Vector)                     :: type_integer_vector
!
  CALL drum !    Open the input and output files
  CALL iosys ('read character options from rwf',-1,0,0,ops)  ! Read in the options string and the Variables
  iosys_on = logkey(ops,'iosys_on',.false.,' ')
  IF ( dollar('$angular_quadrature',card,cpass,inp) ) THEN
       write(iout,1)
       write(iout,2)
       write(iout,1)
!
       type_quadrature = chrkey(card,'type_quadrature','lebedev',' ')
       IF ( type_quadrature == 'lebedev') THEN
            t_start = secnds(0.0)
            test_all = logkey(card,'test_all',.false.,' ')
            print_leb(1)=logkey(card,'print=lebedev=xyz_and_weights',.false.,' ')
            print_leb(2)=logkey(card,'print=lebedev=angles',.false.,' ')
            Call Print_Matrix(type_integer_vector,Rule_Order_Table,title='Lebedev Rules')      
            IF (test_all == .true. ) THEN
                number_of_rules_to_test = rule_max
                DO i = 1, number_of_rules_to_test
!                  Test if available
                   order = rule_order_table(i)
                   avail=.false.
                   IF ( rule_logic_table(i) == 1 ) THEN
                        avail = .true.
                        precision = precision_table (i)
                        write (iout,3) rule_order_table(i), precision
                        Call  Generate_Lebedev_Points_Weights(leb_ang, leb_rule)
                        Call Test_Rule(leb_ang,i)
                        DEALLOCATE(leb_ang%q,leb_ang%w)
                   ELSE              
                        write(iout,4) order
                   END IF
               END DO
            ELSE
               number_of_rules_to_test = intkey(card,'number_of_rules_to_test',1,' ')
               print_leb(1)=logkey(card,'print=lebedev=xyz_and_weights',.false.,' ')
               print_leb(2)=logkey(card,'print=lebedev=angles',.false.,' ')
               ALLOCATE(rule_list(1:number_of_rules_to_test),rule_test(1:rule_max))
               Call intarr(card,'rules',rule_list,number_of_rules_to_test,' ')
               DO i = 1, number_of_rules_to_test
!                         Test if available
                  rule_test(:)= abs(rule_order_table(:) - rule_list(i))
                  iloc = minloc(rule_test,1)
                  order = rule_order_table(iloc)
                  precision = precision_table (iloc)
                  write (iout,3) order, precision
                  Call  Generate_Lebedev_Points_Weights(leb_ang, leb_rule)
                  Call Test_Rule(leb_ang,i)
                  DEALLOCATE(leb_ang%q,leb_ang%w)
               END DO
               DEALLOCATE(rule_list,rule_test)
            END IF
            t_end = secnds(0.0)-t_start       
            write(iout,*) 'Time to compute lebedev angular points and weights = ',t_end
       ELSE IF ( type_quadrature == 'gauss') THEN
            t_start=secnds(0.0)
            theta_ang%n_pts = intkey(card,'number_of_theta_points',10,' ')
            phi_ang%n_pts = intkey(card,'number_of_phi_points',9,' ')
            theta_ang%type_quadrature = chrkey(card,'theta_quadrature_type','gauss',' ')
            phi_ang%type_quadrature = chrkey(card,'phi_quadrature_type','gauss',' ')
            theta_ang%fixed_point=0
            phi_ang%fixed_point=0
            print_gauss(1)=logkey(card,'print=gauss=points_and_weights',.false.,' ')
            Call Gauss_Angular_Grid (theta_ang)
            Call Gauss_Angular_Grid (phi_ang)
            t_end = secnds(0.0) - t_start
            write(iout,*) 'Time to compute Gauss angular points and weights = ',t_end
       ELSE
            Call lnkerr('Quadrature not available')
       END IF
  ELSE IF ( dollar('$atomic_grids',card,cpass,inp) ) THEN
           t_start = secnds(0.0)
           write(iout,1)
           write(iout,8)
           write(iout,1)
           CALL iosys ('read character "grid filename" from rwf',-1,0,0,grid_filename) ! Open the grid file
           CALL iosys ('open grid as new',0,0,0,grid_filename)
! For each center or atom lay out a grid.  
! The grid consists of shells.  In each shell there is an angular and radial set of points.  
! The distribution of both the angular and radial grid is arbitrary.
! The angular grid can be either a non-separable Lebedev grid or a product of theta and phi points.
! If the grid is a lebedev grid then the distribution of points cannot be arbitrary.  
! For a product grids in theta and phi, the breakup into nregions is arbitrary.
! Note that these coordinates are defined with respect to the atomic centers.  We finally need
! them expressed in a single coordinate system.  That is done at a later stage.
           fedvr_functions = logkey(ops,'fedvr_functions',.false.,' ')
           ang_grid = chrkey(ops,'angular_grid','yes',' ')
           rad_grid = chrkey(ops,'radial_grid','yes',' ')
           print_leb(1)=logkey(card,'print=lebedev=xyz_and_weights',.false.,' ')
           print_leb(2)=logkey(card,'print=lebedev=angles',.false.,' ')
           print_ang(1)=logkey(card,'print=angular=points_and_weights',.false.,' ')
           print_ang(2)=logkey(card,'print=angular=polynomials',.false.,' ')
           print_ang(3)=logkey(card,'print=angular=factors',.false.,' ')
           print_rad(1)=logkey(card,'print=radial=points_and_weights',.false.,' ')
           print_rad(2)=logkey(card,'print=radial=polynomials',.false.,' ')
           print_rad(3)=logkey(card,'print=radial=factors',.false.,' ')
           print_norm(1)=logkey(card,'print=normalization=norms',.false.,' ')
           print_norm(2)=logkey(card,'print=normalization=polynomials',.false.,' ')
           print_norm(3)=logkey(ops,'print=normalization=sectors',.false.,' ')
           print_gauss(1)=logkey(card,'print=gauss=points_and_weights',.false.,' ')
           print_gauss(2)=logkey(card,'print=gauss=angles',.false.,' ')
           niter=intkey(ops,'number_of_voronoi_iterations',3,' ')
           alpha=fpkey(ops,'exponential_scattering_grid_cutoff', 1.d0,' ')
           no_voroni=logkey(ops,'voronoi=off',.false.,' ')
           cutoff=fpkey(ops,'scattering_grid_cutoff',10.d0,' ')
           no_scat=logkey(ops,'no_scattering_center',.false.,' ')  ! Default is scattering center present
           one_grid=logkey(ops,'only_one_grid',.false.,' ')        ! Default is multiple grids
           make_3d_grid=logkey(ops,'make_3d_grid',.false.,' ')     ! Default is no
           IF (one_grid) THEN  ! If only one grid, that grid is the scattering grid centered at (0,.0.,0.)
               no_scat=.false.
               no_voroni=.true. ! No Voronoi transformation
           END IF
           no_disk=logkey(ops,'no_disk_output',.false.,' ')
           yukawa_on=logkey(ops,'yukawa=on',.false.,' ')
           yn='yes'
           IF (no_scat==.true.) THEN
               yn='true'
           END IF
           ncent=intkey(card,'number_of_atomic_centers',1,' ')
           ncplus=ncent
           IF (no_scat==.false.) THEN
               ncplus = ncent + 1
           END IF            
           Write(iout,9) no_voroni, niter, no_scat, alpha, cutoff, one_grid, no_disk, yukawa_on 
           IF (iosys_on == .true. ) THEN
               CALL iosys ('write integer "number of atomic centers" to grid',1,ncplus,0,' ')
               CALL iosys ('create real "yukawa exponents" on grid',ncplus,0,0,' ')
           END IF
           ALLOCATE(atom(ncplus),center(ncplus))
           DO i = 1, ncplus
              str='center_'//itoc(i)
              call pakstr(str,len)
              Call fparr(card,'center-'//str(1:len),center(i)%cen,3,' ')
              center(i)%znuc = fpkey(card,'charge-'//str(1:len),1.d0,' ')
              center(i)%R_max = fpkey(card,'atomic_size-'//str(1:len),1.d0,' ')
              center(i)%eta = fpkey(card,'yukawa_exponent-'//str(1:len),1.d0,' ')
           END DO
           t_end = secnds(0.0) - t_start
           write(iout,*) 'Time to set up for atomic grid generation = ',t_end
           DO i=1,ncplus
              str='center_'//itoc(i)
              call pakstr(str,len)
              IF ( dollar('$'//str(1:len),card,cpass,inp) ) THEN
                  atom(i)%n_shl = intkey(card,'number_of_shells',1,' ')
                  atom(i)%znuc = center(i)%znuc
                  atom(i)%R_max = center(i)%R_max
                  atom(i)%cen(:) = center(i)%cen(:)
                  fixed_angular_quadrature=logkey(card,'fixed_angular_quadrature',.false.,' ')
                  IF (iosys_on == .true. ) THEN
                      CALL iosys ('write real "atomic center positions'//str//'" to grid',3,atom(i)%cen,0,' ')
                      CALL iosys ('write real "atomic charge'//str//'" to grid',1,atom%znuc,0,' ')
                  END IF
                  ALLOCATE(atom(i)%shl(1:atom(i)%n_shl))
                  IF (type ==.true.) THEN
                       Write(iout,*)
                       Write(iout,*) '                         Scattering  Center = ', i 
                       Write(iout,*)
                       write(iout,*) '                         Number of Shells   = ', atom(i)%n_shl
                       Write(iout,*)
                       write(iout,*) '                         Fixed_Angular_Grid   = ', fixed_angular_quadrature
                       write(iout,10) atom(i)%cen
                       Call Satshl(atom(i),center)
                  ELSE
                       Write(iout,*)
                       Write(iout,*) '                         Atomic Center      = ', i
                       Write(iout,*)
                       write(iout,*) '                         Number of Shells   = ', atom(i)%n_shl
                       Write(iout,*)
                       write(iout,*) '                         Fixed_Angular_Grid   = ', fixed_angular_quadrature
                       write(iout,11) atom(i)%cen, atom(i)%znuc
                       Call Satshl(atom(i),center) 
                  END IF
              END IF
           END DO
  ELSE
       Call lnkerr('No Input')
  END IF
  Call chainx(0)
1 FORMAT('**********************************************************************************************')
2 FORMAT(/,30x,'               Lebedev/Gauss Quadrature Program               ',/)
3 FORMAT(/,10X,'*** Lebedev Rule  of Order = ',i5,' and Precision = ',i5,' Exists and will be Tested ***')
4 FORMAT(/,10X,'*** Lebedev Rule  of Order = ',i5,' is not available in this program ***')
5 FORMAT(/,10X,'*** Quadrature Rule for Theta Quadrature = ',A16,'  Order = ', i5,' *****', &
        /, 10X,'*** Quadrature Rule for Phi Quadrature   = ',A16,'  Order = ',i5, ' *****')
6 Format(/,25x,' Quadrature Data for Atom = ',i2)
7 Format(/,25x,' Quadrature Data for Scattering Center')
8 FORMAT(/,30x,'               Compute Atomic Grids               ',/)
9 FORMAT(/,10X,'*** Voronoi Off or On = ',l1,'      Number of Voronoi Iterations = ',i5,   &
               '    Scattering Center Present = ',l1,' ***',                               &
         /,10x,'*** Exponential Scattering Grid Cutoff  = ',e15.8,                         &
               '        Scattering Grid Cutoff  = ',e15.8,' ***',                          &
         /,10x,'*** Only One Grid = ',l1,' Disk Output = ',l1,'  Yukawa on = ',l1,' ***')
10 FORMAT(/,10x,'*** Scattering Center Coordinates ***  = ',3(1x,e15.8,1x))
11 FORMAT(/,10x,'*** Atomic Center Coordinates ***  = ',3(1x,e15.8,1x),                    &
          /,10x,'*** Center Charge ***  = ',e15.8)

  Call exit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END PROGRAM m6203
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
