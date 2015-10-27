!deck   m6205
!**begin prologue     m6205
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
  USE Shell_Info
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
  IF ( dollar('$angular_quadrature',card,cpass,inp) ) THEN
       write(iout,1)
       write(iout,2)
       write(iout,1)
!
       type_quadrature = chrkey(card,'type_quadrature','lebedev',' ')
       IF ( type_quadrature == 'lebedev') THEN
            test_all = logkey(card,'test_all',.false.,' ')
            print(1)=logkey(card,'print=XYZ_and_Weights',.false.,' ')
            print(2)=logkey(card,'print=Angles',.false.,' ')
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
               print(1)=logkey(card,'print=XYZ_and_Weights',.false.,' ')
                  print(2)=logkey(card,'print=Angles',.false.,' ')
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
       ELSE IF ( type_quadrature == 'gauss') THEN
            theta_ang%n_pts = intkey(card,'number_of_theta_points',10,' ')
            phi_ang%n_pts = intkey(card,'number_of_phi_points',9,' ')
            theta_ang%type_quadrature = chrkey(card,'theta_quadrature_type','gauss',' ')
            phi_ang%type_quadrature = chrkey(card,'phi_quadrature_type','gauss',' ')
            theta_ang%fixed_point=0
            phi_ang%fixed_point=0
            print(1)=logkey(card,'print=Angles',.false.,' ')
            Call Gauss_Angular_Grid (theta_ang)
            Call Gauss_Angular_Grid (phi_ang)
       ELSE
            Call lnkerr('Quadrature not available')
       END IF
  ELSE IF ( dollar('$atomic_grids',card,cpass,inp) ) THEN
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
! For a product gris intheta and phi, the breakup into nregions is arbitrary.
           only_radial_grid=logkey(ops,'only_radial_grid',.false.,' ')
           only_angular_grid=logkey(ops,'only_angular_grid',.false.,' ')
           niter=intkey(ops,'number_of_voronoi_iterations',3,' ')
           alpha=fpkey(ops,'exponential_scattering_grid_cutoff', 1.d0,' ')
           print(1)=logkey(ops,'print',.false.,' ')
           prnt=.true.
           no_voroni=logkey(ops,'voronoi=off',.false.,' ')
           cutoff=fpkey(ops,'scattering_grid_cutoff',10.d0,' ')
           no_scat=logkey(ops,'no_scattering_center',.false.,' ')  ! Default is scattering center present
           one_grid=logkey(ops,'only_one_grid',.false.,' ')         ! Default is multiple grids
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
           CALL iosys ('write integer "number of atomic centers" to grid',1,ncent,0,' ')
           CALL iosys ('create real "atomic center positions" on grid',3*ncent,0,0,' ')
           CALL iosys ('create real "yukawa exponents" on grid',ncent,0,0,' ')
           CALL iosys ('create real "nuclear charges" on grid',ncent,0,0,' ')
           ALLOCATE(atom(ncplus))
           DO i=1,ncplus
              str='center_'//itoc(i)
              call pakstr(str,len)
              IF ( dollar('$'//str(1:len),card,cpass,inp) ) THEN
                  atom(i)%n_shl = intkey(card,'number_of_shells',1,' ')
                  ALLOCATE(atom(i)%shl(1:atom(i)%n_shl))
                  IF (type ==.true.) THEN
                       Write(iout,*)
                       Write(iout,*) '                         Scattering  Center = ', i 
                       Write(iout,*)
                       write(iout,*) '                         Number of Shells   = ', atom(i)%n_shl
                       Call Satshl(atom(i))
                  ELSE
                       Write(iout,*)
                       Write(iout,*) '                         Atomic Center      = ', i
                       Write(iout,*)
                       write(iout,*) '                         Number of Shells   = ', atom(i)%n_shl
                       Call Satshl(atom(i))
       
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

  stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END PROGRAM m6203
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
