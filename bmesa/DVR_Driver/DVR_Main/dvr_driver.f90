!deck dvr_driver
!**begin prologue     dvr_driver
!**date written       060711   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           dvr!**author             schneider, barry (nsf)
!**source
!**purpose            main dvr driver
!**description        calls routines to input data, construct dvr matrices
!**                   
!**references
!**routines called      Name                    Location    
!                       ----                    --------
!                     Input_Prop             Prop_Sources
!                     Space_Prop             Prop_Sources
!                     Regional_Matrices      Prop_Modules:regional_module
!
!**modules used         Name                    Location               
!                       ----                    --------
!                    dvrprop_global          Modules   
!                    Regional_Matrices       Prop_Modules
!
!**end prologue       dvr_driver
  PROGRAM dvr_driver
  USE dvrprop_global
  USE regional_module
  IMPLICIT NONE
  REAL*4                                   :: secnds
  REAL*4, DIMENSION(10)                    :: del_t
  INTEGER                                  :: iostat 
  INTEGER                                  :: i 
  INTEGER                                  :: j 
  INTEGER                                  :: k 
  INTEGER                                  :: input, output
  INTEGER                                  :: Number_of_Data_Sets
  INTEGER                                  :: intkey
  LOGICAL                                  :: dollar
  LOGICAL                                  :: logkey
  LOGICAL                                  :: Compute_Regional_Matrices
  COMMON /io/ input, output
!
!  Get the input and output file numbers which appear in the 
!  Module input_output.f90 and put them into input and output which appears
!  in the ONLY common block in the code.  This is needed to pass into library
!  routines and for no other purpose.
!
  input = inp
  output = iout
!
!  Open the input and output files
!
  OPEN(input,file='dvr.inp',status='old')
  OPEN(output,file='dvr.out',status='unknown')
!
  IF ( dollar('$begin_data',card,cpass,inp) ) then
       Number_of_Data_Sets = intkey(card,'number_of_data_sets',1,' ')
       Compute_Regional_Matrices = logkey(card,'compute_regional_matrices',.false.,' ')
  END IF
  DO i = 1, number_of_data_sets
     write(iout,1) i
     time(1)=secnds(0.0)
     CALL Input_DVR(i)
     time(2)=secnds(0.0)
     del_t(1) = time(2) - time(1)
!
!    Calculate the Global FEDVR Matrices
!
     CALL Space_DVR(i)
     time(3)=secnds(0.0)
     del_t(2) = time(3) - time(2)
!
!    Compute the Regional Matrices from the FEDVR Global Matrices
!
     IF (Compute_Regional_Matrices ) THEN
         CALL Regional_Matrices
         DO j=1,spdim
            DO k=1,n_reg_real(j)
               DEALLOCATE( mat_reg(k,j)%pt_d, mat_reg(k,j)%ke_mat_d )
            END DO
         END DO
         DEALLOCATE(mat_reg)
         time(4)=secnds(0.0)
         del_t(3) = time(4)-time(3)
     END IF
     DO j=1,spdim
        IF(diag) then                      
           DEALLOCATE(grid(j)%eigv_0,                        &   
                      grid(j)%eigvec_0,                      &   
                      grid(j)%eigv,                          &   
                      grid(j)%eigvec)                    
        END IF
        DEALLOCATE(grid(j)%pt,                               &
                   grid(j)%wt,                               &
                   grid(j)%f,                                &
                   grid(j)%ke,                               &   
                   grid(j)%v)
        DEALLOCATE(grid, num_reg, nfun_reg)
     END DO
     Write(iout,2) 
     Write(iout,3) del_t(1:3)
     Write(iout,2)
  END DO
  CLOSE(input)
  CLOSE(output)
!
1 FORMAT(/,25x,'Data Set = ',i3)
2 FORMAT('*****************************************************************')
3 FORMAT( /,10X,'Time to input basic data               = ',f15.8,/,     &
          /,10X,'Time to compute global FEDVR spatial ',                 &
          /,10x,'matrices                               = ',f15.8,/,     &
          /,10x,'Time to compute the regional spatial '                  &
          /,10x,'matrices                               = ',f15.8)
  stop
END PROGRAM dvr_driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
