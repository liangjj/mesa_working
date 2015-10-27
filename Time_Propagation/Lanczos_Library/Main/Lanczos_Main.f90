!
!deck Lanczos_Main
!**begin prologue     Lanczos_Main
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**description        
!**                   
!**                   
!**references
!**routines called      Name                    Location    
!                       ----                    --------
!
!**modules used         Name                    Location               
!                       ----                    --------
!
!**end prologue       Lanczos_Main
  PROGRAM Lanczos_Main
  USE dvrprop_global
  USE Lanczos_Driver_Module
  USE Matrix_Module
  USE Lanczos_Module
!
  IMPLICIT NONE
  REAL*4                                   :: secnds
  REAL*4, DIMENSION(10)                    :: del_t
  INTEGER                                  :: iostat 
  CHARACTER (LEN=8)                        :: matrix_source
!
! Read in the Input
!
  time(1)=secnds(0.0)
  CALL Lanczos_Input(matrix_source)
  time(2)=secnds(0.0)
  del_t(1) = time(2) - time(1)
  WRITE(iout,1)
  WRITE(iout,2) del_t(1)
  WRITE(iout,1)
!
!
  CALL Read_Matrix
  time(3)=secnds(0.0)
  del_t(2) = time(3) - time(2)
  WRITE(iout,1)
  WRITE(iout,3) del_t(2)
  WRITE(iout,1)
  IF ( .not.solve_only ) THEN
!
     time(1)=secnds(0.0)
!
     CALL Iterative_Data
!        
!      
     CALL Lanczos_Driver
     time(2) = secnds(0.0)
     del_t(1) = time(2) - time(1)
     WRITE(iout,1)
     WRITE(iout,4) del_t(1)
     WRITE(iout,1)
  END IF
1 FORMAT('***********************************************'                           &
         '*************************')
2 FORMAT(/,10X,'Time to Input Basic Data               = ',f15.8)
3 FORMAT(/,10X,'Time to Read, Factor and Diagonalize Matrices = ',f15.8)
4 FORMAT(/,10x,'Time for the Lanczos Procedure                = ',f15.8)
  stop
END PROGRAM Lanczos_Main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
