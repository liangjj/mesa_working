!deck   Angular_Driver
!**begin prologue     Angular_Driver
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
  PROGRAM Angular_Driver
  USE Lebedev_Quadrature
  USE Gauss_Angular_Quadrature
  IMPLICIT NONE
  CHARACTER (LEN=3)                        :: itoc
  CHARACTER (LEN=80)                       :: chrkey
  LOGICAL                                  :: dollar
  LOGICAL                                  :: logkey
  LOGICAL                                  :: avail
  LOGICAL                                  :: test_all
  REAL(idp)                                :: fpkey
  REAL(idp), DIMENSION(4)                  :: range
  INTEGER                                  :: intkey
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: iloc
  INTEGER                                  :: st
  INTEGER                                  :: irange
  INTEGER, DIMENSION(:), ALLOCATABLE       :: rule_list
  INTEGER, DIMENSION(:), ALLOCATABLE       :: rule_test
  INTEGER                                  :: number_of_rules_to_test
  CHARACTER (LEN = 80)                     :: title
  CHARACTER (LEN = 80)                     :: type_quadrature
!  TYPE(Real_Vector)                        :: type_real_vector
!  TYPE(Real_Matrix)                        :: type_real_matrix
  TYPE(Integer_Vector)                     :: type_integer_vector
  TYPE(LEBEDEV)                            :: leb_ang
  TYPE(THETA)                              :: theta_ang
  TYPE(PHI)                                :: phi_ang
!
  open (inp,file='Ang.input',access='sequential',form='formatted',iostat=st,status='old')
  open (iout,file='Ang.output',access='sequential',form='formatted',iostat=st,status='unknown')
  CALL drum !    Open the input and output files
  CALL iosys ('read character options from rwf',-1,0,0,ops)  ! Read in the options string and the Variables
  write(iout,1)
  write(iout,2)
  write(iout,1)
  IF ( dollar('$angular_test',card,cpass,inp) ) THEN
!
       type_quadrature = chrkey(card,'type_quadrature','lebedev',' ')
       IF ( type_quadrature == 'lebedev') THEN
            test_all = logkey(card,'test_all',.false.,' ')
            print=logkey(card,'print=on',.false.,' ')
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
               print=logkey(card,'print=on',.false.,' ')
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
                  IF (print == .true. ) THEN
                      Call Print_Matrix(type_real_matrix,leb_ang%q,3,order,title='Points',frmt='fr')
                      Call Print_Matrix(type_real_vector,leb_ang%w,title='Weights',frmt='fr')
                  END IF
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
            Call Gauss_Angular_Grid (theta_ang)
            Call Gauss_Angular_Grid (phi_ang)
            print=.true.
            IF (print == .true. ) THEN
                Call Print_Matrix(type_real_vector,theta_ang%q,title='Cos Theta',frmt='fr')
                Call Print_Matrix(type_real_vector,theta_ang%sin_thet,title='Sin Theta',frmt='fr')
                Call Print_Matrix(type_real_vector,theta_ang%wt,title='Theta Weights',frmt='fr')
                Call Print_Matrix(type_real_vector,phi_ang%q,title='Phi',frmt='fr')
                Call Print_Matrix(type_real_vector,phi_ang%cos_phi,title='Cos Phi',frmt='fr')
                Call Print_Matrix(type_real_vector,phi_ang%sin_phi,title='Sin Phi',frmt='fr')
                Call Print_Matrix(type_real_vector,phi_ang%wt,title='Phi Weights',frmt='fr')
            END IF
       ELSE
            Call lnkerr('Quadrature not available')
       END IF
  ELSE
       Call lnkerr('No Input')
  END IF
1 FORMAT('***********************************************'                           &
         '*************************')
2 FORMAT(/,30x,'*** Lebedev Quadrature Test Program ***')
3 FORMAT(/,10X,'*** Lebedev Rule  of Order = ',i5,' and Precision = ',i5,' Exists and will be Tested ***')
4 FORMAT(/,10X,'*** Lebedev Rule  of Order = ',i5,' is not available in this program ***')
5 FORMAT(/,10X,'*** Quadrature Rule for Theta Quadrature = ',A16,'  Order = ', i5,' *****', &
        /, 10X,'*** Quadrature Rule for Phi Quadrature   = ',A16,'  Order = ',i5, ' *****')
  stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END PROGRAM Angular_Driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
