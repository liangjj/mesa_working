!***********************************************************************
! Lanczos_Driver_Module
!**begin prologue     Lanczos_Driver_Module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Arnoldi, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose
!***  
!***description 
!*** 
!**references
!**modules needed     See USE statements below
!**comments           
!**                   
!**                   
!**                   
!**                   
!**end prologue       Lanczos_Driver_Module
!***********************************************************************
!***********************************************************************
                           MODULE Lanczos_Driver_Module
                           USE Matrix_Module
                           USE Lanczos_Module
                           USE Iterative_Global
                           USE initial_state_module
                           USE Pack_Global
                           USE Pack_Matrix_Module
!***********************************************************************
!***********************************************************************
                           INTERFACE drive_lanczos
             MODULE PROCEDURE drive_lanczos_d,                           &
                              drive_lanczos_h
                       END INTERFACE drive_lanczos
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!**********************************************************************
!                            Driving Routine
!**********************************************************************
!***********************************************************************
! lanczos_driver
!**begin prologue     lanczos_driver
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             
!**purpose            
!***                  
!***                  
!***                  
!***description 
!***            
!***            
!***            
!***            
!***            
!***            
!***            
!***            
!**references
!**modules needed     See USE statements below
!**end prologue       lanczos_driver
!***********************************************************************
!***********************************************************************
  SUBROUTINE lanczos_driver
  IMPLICIT NONE
  INTEGER                       :: i 
  INTEGER                       :: IOSTAT 
  lwork = 5*(maxvec+1)
!
! Allocate memory for the main arrays.              
!                                                                         
! Get Data on Initial State
!
  OPEN(UNIT=99,FILE='wavefunction',ACCESS='sequential',        &
           FORM='unformatted',IOSTAT=IOSTAT,STATUS='unknown')
  CALL Initial_State_Data
!
  IF ( matrix_type == 'real' ) THEN
       ALLOCATE( psi_d(1:n3d) )
       ALLOCATE( matrix_diagonal_d(1:n3d), packed_columns_v_d(1:lenbuf),   &
                 non_zero_columns_v(1:n3d), row_index_v(1:lenbuf) )
       Call Pack_Matrix(matrix_d,matrix_diagonal_d,packed_columns_v_d,     &
                        non_zero_columns_v, row_index_v,n3d)
       CALL Initial_Vector(psi_d)
       CALL Drive_Lanczos(psi_d )
       DEALLOCATE( matrix_diagonal_d, packed_columns_v_d, non_zero_columns_v, row_index_v )
       DEALLOCATE( psi_d, rhs_d )
  ELSE IF(matrix_type == 'hermitian' ) THEN
       ALLOCATE( psi_z(1:n3d) )
       ALLOCATE(matrix_diagonal_z(1:n3d), packed_columns_v_z(1:lenbuf),    &
                non_zero_columns_v(1:n3d), row_index_v(1:lenbuf) )
       Call Pack_Matrix(matrix_z,matrix_diagonal_z,packed_columns_v_z,     &
                        non_zero_columns_v,row_index_v,n3d)
       CALL Initial_Vector(psi_z)
       CALL Drive_Lanczos(psi_z )
       DEALLOCATE( matrix_diagonal_z, packed_columns_v_z, non_zero_columns_v, row_index_v )
       DEALLOCATE( psi_z, rhs_z )
  ELSE
       call lnkerr('Quit.  Bad Keyword')
  END IF
  CLOSE(UNIT=99)
END SUBROUTINE lanczos_Driver
!**********************************************************************
!**********************************************************************
!deck Drive_lanczos_d
!**begin prologue     Drive_Lanczos_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords
!**
!**author             schneider, b. i.(nsf)
!**source  
!**purpose 
!**        
!**        
!**        
!**references         Drive_Lanczos_d
!**routines called    iosys, util and mdutil
!**end prologue       
  SUBROUTINE Drive_Lanczos_d(wave_function)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                     :: wave_function
!
! Allocate Memory
!
  WRITE(iout,1)
  WRITE(iout,2)
  ALLOCATE( vec_d(n3d,0:maxvec), h_vec_d(n3d), h_mat_d(0:maxvec,0:2),                                         &
            rhs_tri_d(0:maxvec), soln_tri_d(0:maxvec), lanczos_tri_d(0:maxvec), a(0:maxvec), b(0:maxvec),     &
            work_d(0:lwork) )
  WRITE(iout,1)
  CALL lanczos(wave_function)
!
! Deallocate Memory
  WRITE(iout,3)
!
  DEALLOCATE( vec_d, h_vec_d, h_mat_d, rhs_tri_d, soln_tri_d, lanczos_tri_d, a, b, work_d )
  WRITE(iout,1)
!
1 FORMAT('***********************************************'                                                     &
         '*************************')
2 FORMAT(/,20x,'Begin Lanczos Allocate Storage')
3 FORMAT(/,20x,'Converged: End Allocate Storage')
END SUBROUTINE Drive_Lanczos_d
!***********************************************************************
!***********************************************************************
!deck Drive_Lanczos_h
!**begin prologue     Drive_Lanczos_h
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords         
!**
!**author             schneider, b. i.(nsf)
!**source           
!**purpose          
!**                 
!**                 
!**                 
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Drive_Lanczos_h
  SUBROUTINE Drive_Lanczos_h(wave_function)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)                 :: wave_function                            
!
! Allocate Memory
!
  WRITE(iout,1)
  WRITE(iout,2)
  ALLOCATE( vec_z(n3d,0:maxvec), h_vec_z(n3d), h_mat_z(0:maxvec,0:maxvec),                                    &
            rhs_tri_z(0:maxvec), soln_tri_z(0:maxvec), lanczos_tri_z(0:maxvec), a(0:maxvec), b(0:maxvec),     &
            work_z(0:lwork) )
  WRITE(iout,1)
  CALL lanczos(wave_function)
!
! Deallocate Memory
  WRITE(iout,3)
!
  DEALLOCATE( vec_z, h_vec_z, h_mat_z, rhs_tri_z, soln_tri_z, lanczos_tri_z, a, b, work_z )
  WRITE(iout,1)
!
1 FORMAT('***********************************************'                                   &
         '*************************')
2 FORMAT(/,20x,'Begin Lanczos Allocate Storage')
3 FORMAT(/,20x,'Converged: End Allocate Storage')
  WRITE(iout,1)
END SUBROUTINE Drive_Lanczos_h
!***********************************************************************
!***********************************************************************
END MODULE Lanczos_Driver_Module
!***********************************************************************
!***********************************************************************
