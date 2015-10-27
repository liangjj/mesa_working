!
  MODULE Time_Propagation_Module
  USE Data
  USE Matrix_Print
  IMPLICIT NONE
!
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                            INTERFACE Time_Propagation                       
                       MODULE PROCEDURE Propagation_CN_Length,            &
                                        Propagation_CN_Velocity,          &
                                        Real_Propagation_SO_Length,       &
                                        Complex_Propagation_SO_Length
                            END INTERFACE Time_Propagation
!
                            INTERFACE Propagators                       
                       MODULE PROCEDURE Real_Propagator,                  &
                                        Complex_Propagator  
                            END INTERFACE Propagators
!____________________________________________________________________________________________!
!____________________________________________________________________________________________!

  TYPE MATRICES
     COMPLEX(idp), DIMENSION(:,:), ALLOCATABLE         :: exp_real
     REAL(idp), DIMENSION(:,:), ALLOCATABLE            :: exp_imag
  END TYPE MATRICES

  TYPE REAL_TIME
     TYPE(MATRICES), DIMENSION(:), ALLOCATABLE         :: mat       
  END TYPE REAL_TIME

  TYPE (REAL_TIME)                                     :: r_prop

  TYPE IMAG_TIME
     TYPE(MATRICES), DIMENSION(:), ALLOCATABLE         :: mat       
  END TYPE IMAG_TIME

  TYPE (IMAG_TIME)                                     :: i_prop
!
  TYPE CN_LENGTH
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: ham_l
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: ham_u
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: ham_d
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: sol
  END TYPE CN_LENGTH
!
!
  TYPE CN_VELOCITY
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: ham_l
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: ham_u
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: ham_d
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: sol
  END TYPE CN_VELOCITY
!
  TYPE REAL_SO_LENGTH
     REAL(idp), DIMENSION(:), ALLOCATABLE              :: ham_l
     REAL(idp), DIMENSION(:), ALLOCATABLE              :: ham_u
     REAL(idp), DIMENSION(:), ALLOCATABLE              :: ham_d
     REAL(idp), DIMENSION(:), ALLOCATABLE              :: sol
  END TYPE REAL_SO_LENGTH
!
  TYPE COMPLEX_SO_LENGTH
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: ham_l
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: ham_u
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: ham_d
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: sol
  END TYPE COMPLEX_SO_LENGTH
!
  TYPE(CN_LENGTH)                                      :: cn_len
!
  TYPE(CN_VELOCITY)                                    :: cn_vel
!
  TYPE(REAL_SO_LENGTH)                                 :: real_so_len
!
  TYPE(COMPLEX_SO_LENGTH)                              :: complex_so_len

                                       Contains
!********************************************************************************
!********************************************************************************
  SUBROUTINE Real_Propagation_SO_Length (real_so_len)
  IMPLICIT NONE
  TYPE(REAL_SO_LENGTH)                       :: real_so_len
  INTEGER                                    :: i
  INTEGER                                    :: INFO
!********************************************************************************
  END SUBROUTINE Real_Propagation_SO_Length
!********************************************************************************
!********************************************************************************
!********************************************************************************
  SUBROUTINE Complex_Propagation_SO_Length (complex_so_len)
  IMPLICIT NONE
  TYPE(COMPLEX_SO_LENGTH)                       :: complex_so_len
  INTEGER                                    :: i
  INTEGER                                    :: INFO
!********************************************************************************
  END SUBROUTINE Real_Propagation_SO_Length
!********************************************************************************
  SUBROUTINE COMPLEX_PROPAGATOR(r_prop)
  IMPLICIT NONE
  TYPE (REAL_TIME)                               :: r_prop
  TYPE(MATRICES), DIMENSION(:), ALLOCATABLE      :: mat       
  INTEGER                                        :: n
  INTEGER                                        :: i
  INTEGER                                        :: j
  INTEGER                                        :: k
  COMPLEX(idp), DIMENSION(2), ALLOCATABLE        :: EigVal
  COMPLEX(idp), DIMENSION(2,2), ALLOCATABLE      :: EigVec
  CHARACTER (LEN=5)                              :: itoc 
!
  ALLOCATE(r_prop%mat(1:M_Size))
  DO i =1, M_Size
     ALLOCATE(r_prop%mat(i)%exp_real(2,2))
  END DO
!
  DO i = 2, M_Size-1
     EigVec(1,1) = D(i)*.5d0
     EigVec(2,2) = D(i)*.5d0
     EigVec(1,2) = E(i)
     EigVec(2,1) = E(i)


                              ! FIRST BLOCK !
                     ! FIST 2X2 BLOCK IS AN ODD BLOCK !

  CALL eigV(D(1),E(1),E(1),D(2)*0.5d0,eigVal,eigVec)
  DO i =1, 2
     DO j =1, i
        DO k =1, 2
           Mtrx(i,j,1) = Mtrx(i,j,1) + eigVec(i,k)*exp(-eye*eigVal(k)*delta_t*0.5d0)*eigVec(j,k)
        END DO
     END DO
  END DO

                               ! LAST BLOCK !
             ! IF M_SIZE IS EVEN, LAST 2x2 BLOCK IS AN ODD BLOCK !
             ! IF M_SIZE IS ODD, LAST 2X2 BLOCK IS AN EVEN BLOCK !

  CALL eigV(D(M_Size-1)*0.5d0,E(1),E(1),D(M_Size)*0.5d0,eigVal,eigVec)
  DO i =1, 2
     DO j =1, i
        DO k =1, 2
           Mtrx(i,j,M_Size-1) = Mtrx(i,j,M_Size-1) +eigVec(i,k)*exp(-eye*eigVal(k)*delta_t)*eigVec(j,k)
        END DO
     END DO
  END DO

                              ! MIDDLE BLOCKS !
  DO n = 2, M_Size - 2, 2
                               ! EVEN BLOCKS !
     CALL eigV(D(n)*0.5d0, E(1), E(1),D(n+1)*0.5d0,eigVal,eigVec)
     DO i =1, 2
         DO j =1, i
            DO k =1, 2
               Mtrx(i,j,n) = Mtrx(i,j,n) +eigVec(i,k)*exp(-eye*eigVal(k)*delta_t)*eigVec(j,k)
            END DO
         END DO
      END DO
                               ! ODD BLOCKS !
     CALL eigV(D(n+1)*0.5d0, E(1), E(1), D(n+2)*0.5d0,eigVal,eigVec)
     DO i =1, 2
         DO j =1, i
            DO k =1, 2
               Mtrx(i,j,n+1) = Mtrx(i,j,n+1) + eigVec(i,k)*exp(-eye*eigVal(k)*delta_t*0.5d0)*eigVec(j,k)
            END DO
         END DO
      END DO
  END DO
  Mtrx(1,2,:) = Mtrx(2,1,:)
  
  DEALLOCATE (eigVal, eigVec)
!*******************************************************************************
  END SUBROUTINE COMPLEX_PROPAGATOR
!*******************************************************************************
  SUBROUTINE REAL_PROPAGATOR(i_prop)
  IMPLICIT NONE
  TYPE (IMAG_TIME)                               :: i_prop
  TYPE(MATRICES), DIMENSION(:), ALLOCATABLE      :: mat       
  INTEGER                                        :: n
  INTEGER                                        :: i
  INTEGER                                        :: j
  INTEGER                                        :: k
  REAL(idp), DIMENSION(:), ALLOCATABLE           :: eigVal
  REAL(idp), DIMENSION(:,:), ALLOCATABLE         :: eigVec
!  REAL(idp), DIMENSION(2,2, M_Size-1)            :: Mtrx
  REAL(idp), DIMENSION(:,:,:),ALLOCATABLE        :: Mtrx
  CHARACTER(LEN=5)                               :: itoc 
  TYPE(Real_Vector)                              :: type_real_vector
  TYPE(Real_Matrix)                              :: type_real_matrix

  ALLOCATE (eigVal(2), eigVec(2,2) )

  Mtrx(:,:,:)        = 0.0d0
  
                           ! FIRST BLOCK !
                     ! FIST 2X2 BLOCK IS AN ODD BLOCK !

  CALL eigV(D(1),E(1),E(1),D(2)*0.5d0,eigVal,eigVec)
  DO i =1, 2
     DO j =1, i
        DO k =1, 2
           Mtrx(i,j,1) = Mtrx(i,j,1) + eigVec(i,k)*exp(-eigVal(k)*delta_t*0.5d0)*eigVec(j,k)
        END DO
     END DO
  END DO

                               ! LAST BLOCK !
             ! IF M_SIZE IS EVEN, LAST 2x2 BLOCK IS AN ODD BLOCK !
             ! IF M_SIZE IS ODD, LAST 2X2 BLOCK IS AN EVEN BLOCK !

  CALL eigV(D(M_Size-1)*0.5d0,E(1),E(1),D(M_Size),eigVal,eigVec)
  DO i =1, 2
     DO j =1, i
        DO k =1, 2
           Mtrx(i,j,M_Size-1) = Mtrx(i,j,M_Size-1) +eigVec(i,k)*exp(-eigVal(k)*delta_t)*eigVec(j,k)
        END DO
     END DO
  END DO


                              ! MIDDLE BLOCKS !
  DO n = 2, M_Size - 2, 2
     CALL eigV(D(n)*0.5d0, E(1), E(1),D(n+1)*0.5d0,eigVal,eigVec)
     DO i =1, 2
         DO j =1, i
            DO k =1, 2
               Mtrx(i,j,n) = Mtrx(i,j,n) +eigVec(i,k)*exp(-eigVal(k)*delta_t)*eigVec(j,k)
            END DO
         END DO
     END DO
                               ! ODD BLOCKS !
     CALL eigV(D(n+1)*0.5d0, E(1), E(1), D(n+2)*0.5d0,eigVal,eigVec)
     DO i =1, 2
         DO j =1, i
            DO k =1, 2
               Mtrx(i,j,n+1) = Mtrx(i,j,n+1) + eigVec(i,k)*exp(-eigVal(k)*delta_t*0.5d0)*eigVec(j,k)
            END DO
         END DO
      END DO
  END DO
  Mtrx(1,2,:) = Mtrx(2,1,:)


!  DO i =1, M_Size-1
!     eigVec(:,:) = Mtrx(:,:,i)
!     Call Print_Matrix(type_real_matrix,eigVec,2, 2,title='exponential block ' //itoc(i))  
!  END DO


  DEALLOCATE (eigVal, eigVec)
!*******************************************************************************
  END SUBROUTINE REAL_PROPAGATOR
!*******************************************************************************
  END MODULE Time_Propagation_Module
!********************************************************************************
