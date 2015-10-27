!deck Driver
!**begin prologue     Driver
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**description        
!**                   
!**                   
!**references
!**routines called      Name                    Location    
!                       ----                    --------
!                     Input_Prop             CC_Prop_Sources
!                     Space_Prop             CC_Prop_Sources
!                     Regional_Matrices      CC_Prop_Modules:regional_module
!                     v_couple               CC_Prop_Sources
!                     Propagation_Driver     CC_Prop_Modules:Propagation_Module

!
!**modules used         Name                    Location               
!                       ----                    --------
!                    dvrprop_global          Modules   
!                    Propagation_Module      CC_Prop_Modules
!
!**end prologue       Driver
  PROGRAM Driver
  USE input_output
  USE Davidson_Module
  IMPLICIT NONE
  REAL*4                                   :: secnds
  REAL*4, DIMENSION(10)                    :: del_t
  INTEGER                                  :: iostat 
  INTEGER                                  :: input
  INTEGER                                  :: output
  COMMON /io/ input, output
!
! Read in the Input
!
  input = inp
  output = iout
!
!  Open the input and output files
!
  OPEN(input,file='Davidson.inp',status='old')
  OPEN(output,file='Davidson.out',status='unknown')
  c_time(1)=secnds(0.0)
  CALL Davidson_Data
  c_time(2)=secnds(0.0)
  del_t(1) = c_time(2) - c_time(1)
  WRITE(iout,1)
  WRITE(iout,2) del_t(1)
  WRITE(iout,1)
!
                      Matrix_Input_Section:              &
!
  IF(matrix_source == 'from_input') THEN
!
!
     CALL Read_Matrix_From_Input
     c_time(3)=secnds(0.0)
     del_t(2) = c_time(3) - c_time(2)
!
!    ALLOCATE the Memory for the Matrices used to Store the Wavefunction
!    and for the Potential.  Both the Lanczos and Split Operator need these
!    arrays.
!
  ELSE IF(matrix_source == 'from_disk') THEN 
!    CALL Read_Matrix_From_Disk
     c_time(3)=secnds(0.0)
  END IF              Matrix_Input_Section  
!
  WRITE(iout,1)
  WRITE(iout,3) del_t(2)
  WRITE(iout,1)
!
! Calculate either the guess vectors or the guess right hand side
!
  IF (type_calculation == 'linear_system') THEN
      IF (matrix_type == 'real_symmetric') THEN 
          ALLOCATE(b_d(matrix_size))
          Call guess_solution(eigenvectors_d,b_d)
      ELSE IF (matrix_type == 'hermitian') THEN      
          ALLOCATE(b_z(matrix_size))
          Call guess_solution(eigenvectors_z,b_z)
      END IF
  ELSE IF (type_calculation == 'eigenvalues') THEN
      IF (matrix_type == 'real_symmetric') THEN 
          ALLOCATE(eigenvectors_d(matrix_size,guess_size))
          Call guess_solution(eigenvectors_d,b_d)
      ELSE IF (matrix_type == 'hermitian') THEN      
          ALLOCATE(eigenvectors_z(matrix_size,guess_size))
          Call guess_solution(eigenvectors_z,b_z)
      END IF
  END IF
!  stop
  WRITE(iout,1)
  Call Allocation
  WRITE(iout,1)
  Write(iout,4)
  total_time = 0.0
  ALLOCATE(total_time(20))
  IF (type_calculation == 'linear_system') THEN
      IF (matrix_type == 'real_symmetric') THEN 
          v_in_d = b_d 
          Call linear_system_driver(v_in_d,v_out_d,rhs_d)
      ELSE IF (matrix_type == 'hermitian') THEN      
          v_in_z = b_z
          Call linear_system_driver(v_in_z,v_out_z,rhs_z)
      END IF
  ELSE IF (type_calculation == 'eigenvalues') THEN      
      IF (matrix_type == 'real_symmetric') THEN 
!          Call eigen_system_driver(v_in_d,v_out_d,rhs_d)
      ELSE IF (matrix_type == 'hermitian') THEN      
!          Call eigen_system_driver(v_in_z,v_out_z,rhs_z)
      END IF
  END IF
1 FORMAT('***********************************************'                           &
         '*************************')
2 FORMAT(/,10X,'Time to Input Basic Data               = ',f15.8)
3 FORMAT(/,10X,'Time to Input the Matrices             = ',f15.8)
4 FORMAT(/,10X,'Starting the Main Calculation')
  stop
END PROGRAM Driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
