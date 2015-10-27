!***********************************************************************
! Exponential_on_Vector_Module
!**begin prologue     Exponential_on_Vector_Module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            To perform the operation of the exponentiation of a matrix on a vector e
!***                  when the matrix exists in diagonal form.
!***description       
!***                  
!***                  
!***                  
!***                  
!***                  
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
!**comments           
!**                   
!**                   
!**                   
!**                   
!**end prologue       Exponential_on_Vector_Module
!***********************************************************************
!***********************************************************************
                           MODULE Exponential_on_Vector_Module
                           USE dvr_matrix_vector_multiply_module
                           USE Matrix_Print
                             IMPLICIT NONE
                           TYPE(REAL_MATRIX)      :: type_real_matrix
                           TYPE(COMPLEX_MATRIX)   :: type_complex_matrix
                           TYPE(REAL_VECTOR)      :: type_real_vector
                           TYPE(COMPLEX_VECTOR)   :: type_complex_vector
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                           INTERFACE Exponential_on_Vector
             MODULE PROCEDURE Exponential_on_Vector_d,                            &
                              Exponential_on_Vector_z                             
                       END INTERFACE Exponential_on_Vector
!
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!**begin prologue     Exponential_on_Vector_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose      
!***            
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Exponential_on_Vector_d
!***********************************************************************
  SUBROUTINE Exponential_on_Vector_d(vector_out,vector_in,lanczos_eigenvalues,lanczos_eigenvectors,it)
  IMPLICIT NONE
  REAL(idp), DIMENSION(:)                   :: vector_in
  REAL(idp), DIMENSION(:)                   :: vector_out
  REAL(idp), DIMENSION(0:it)                :: lanczos_eigenvalues
  REAL(idp), DIMENSION(0:maxvec,0:it)       :: lanczos_eigenvectors
  INTEGER                                   :: it
  CHARACTER(LEN=4)                          :: itoc
!
!               
!             Compute the projection of the time-dependent state on the initial vector,
!             The initial vector is the first Lanczos vector.  So, this just scales
!             the first component of each eigenvector by an exponential.
!
!             C(t) = Sum_{lambda} Exp(-lambda * deltat ) d_{lambda}(0) C_{lambda}
!                                                     +
!             Compute the d_{lambda}(0) = [C_{lambda}]  S C(0)
!
  IF(log_iterative(10)) THEN
     call Print_Matrix(type_real_vector,vector_in(1:it+1),                                      &
                       title='Projection of Initial Vector on Eigenbasis iteration = '//itoc(it))
  END IF
! 
!           Here is the scaling.
!
  vector_in(1:it+1) = EXP(-lanczos_eigenvalues(0:it) * deltat/hbar) * vector_in(1:it+1)
  IF(log_iterative(10)) THEN
      call Print_Matrix(type_real_vector,vector_in(1:it+1),                                     &
                        title='Projection Multiplied by Exponential iteration = '//itoc(it))
  END IF
!
!  
!             Now express things in the Lanczos basis instead of the eigenfunction basis.
!
  call ebcx(vector_out,maxvec+1,lanczos_eigenvectors,maxvec+1,vector_in,maxvec+1,it+1,it+1,1)
!
END SUBROUTINE Exponential_on_Vector_d
!***********************************************************************
!***********************************************************************
!**begin prologue     Exponential_on_Vector_z
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose      
!***            
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Exponential_on_Vector_z
!***********************************************************************
  SUBROUTINE Exponential_on_Vector_z(vector_out,vector_in,lanczos_eigenvalues,lanczos_eigenvectors,it)
  IMPLICIT NONE
  COMPLEX(idp), DIMENSION(:)                     :: vector_in
  COMPLEX(idp), DIMENSION(:)                     :: vector_out
  REAL(idp), DIMENSION(0:it)                     :: lanczos_eigenvalues
  REAL(idp), DIMENSION(0:maxvec,0:it)            :: lanczos_eigenvectors
  INTEGER                                        :: it
  CHARACTER(LEN=4)                               :: itoc
!
!               
!             Compute the projection of the time-dependent state on the initial vector,
!             The initial vector is the first Lanczos vector.  So, this just scales
!             the first component of each eigenvector by an exponential.
!
!             C(t) = Sum_{lambda} Exp(-lambda * deltat ) d_{lambda}(0) C_{lambda}
!                                                     +
!             Compute the d_{lambda}(0) = [C_{lambda}]  S C(0)
!
  IF(log_iterative(10)) THEN
     call Print_Matrix(type_complex_vector,vector_in(1:it+1),                                   &
                       title='Projection of Initial Vector on Eigenbasis iteration = '//itoc(it))
  END IF
! 
!           Here is the scaling.
!
  vector_in(1:it+1) = EXP(-eye*lanczos_eigenvalues(0:it) * deltat/hbar) * vector_in(1:it+1)
  IF(log_iterative(10)) THEN
     call Print_Matrix(type_complex_vector,vector_in(1:it+1),                                   &
                       title='Projection Multiplied by Exponential iteration = '//itoc(it))
  END IF
!
!  
!             Now express things in the Lanczos basis instead of the eigenfunction basis.
!
  call ebccx(vector_out,maxvec+1,lanczos_eigenvectors,maxvec+1,vector_in,maxvec+1,it+1,it+1,1)
!
END SUBROUTINE Exponential_on_Vector_z
!***********************************************************************
!***********************************************************************
           END MODULE Exponential_on_Vector_Module
!***********************************************************************
!***********************************************************************
