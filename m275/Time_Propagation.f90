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
                            INTERFACE CN_Propagation                       
                       MODULE PROCEDURE CN_Length,            &
                                        CN_Velocity,          &
                            END INTERFACE CN_Propagation
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
!____________________________________________________________________________________________!
!____________________________________________________________________________________________!
!_______________________________________TIME__INDEPENDENT____________________________________!
!______________________(FINDING EIGENVALUES AND EIGENVECTORS USING DSTEVX)___________________!
!____________________________________________________________________________________________!
!
!                                     DSTEVX:
! If Y_N='N', only eigenvalues will be computed; 
! If Y_N='V', eigenvalues and eigenvectors wil be computed.
!
! If Range = 'A', all eigenvalues will be computed.
! If Range = 'V', all eigenvalues in range (VL,VU] will be computed.
! If Range = 'I', all eigenvalues from ILth through IUth will be computed.
! 

  SUBROUTINE H_0 
  IMPLICIT NONE
  REAL(idp), DIMENSION(:), ALLOCATABLE     :: WORK 
  REAL(idp), DIMENSION(:), ALLOCATABLE     :: IWORK 
  REAL(idp), DIMENSION(:), ALLOCATABLE     :: IFAIL 
  REAL(idp)                                :: XX
  REAL(idp)                                :: AD
  REAL(idp)                                :: ABSTOL
  REAL(idp)                                :: DLAMCH
  INTEGER                                  :: i 
  INTEGER                                  :: VL 
  INTEGER                                  :: VU
  INTEGER                                  :: M_FOUND
  INTEGER                                  :: NSPLIT
  INTEGER                                  :: INFO
  INTEGER, DIMENSION(:), ALLOCATABLE       :: ISPLIT
  INTEGER, DIMENSION(:), ALLOCATABLE       :: IBLOCK
  TYPE(Real_Vector)                        :: type_real_vector
  TYPE(Real_Matrix)                        :: type_real_matrix
!
  ALLOCATE(D(1:M_Size), E(1:M_Size), IBLOCK(1:M_Size), ISPLIT(1:M_Size),      &
           WORK(1:5*M_Size), W(1:M_Size), IWORK(1:5*M_Size), IFAIL(1:M_Size) )
  IF (Get_Eigenvectors == .true.) THEN
      Y_N = 'V'
      ALLOCATE(Z(M_Size,Number_of_Eigenvectors))
  END IF
  ABSTOL = 2 * DLAMCH('S')
  AD = - 1.d0 / ( 2.d0 * Step_Size * Step_Size )
  D(1:M_Size) = - 2.d0 * AD
  E(1:M_Size) = AD
  IF ( Add_Potential  == .true. ) THEN
       ALLOCATE(V(1:M_Size))
       XX = left_end
       DO i = 1, M_Size
          XX = XX + Step_Size
          V(i) = - 1.d0 / Sqrt (  ( 1.d0 + XX * XX ) )
          D(i) = D(i) +V(i) 
       END DO
  END IF
  Call cpu_time(t_i)
  Call DSTEBZ('I', 'E', M_Size, VL, VU, IL, IU, ABSTOL, D, E, M_FOUND, NSPLIT, W,    &
               IBLOCK, ISPLIT, WORK, IWORK, INFO)
  Call cpu_time(t_f)
  write(outdat,*) '               ELapsed Time for eigenvalues= ',t_f-t_i
  IF (INFO /= 0 ) THEN
      Write(outdat,*) '           Eigenvalue Routine Failed'
  ELSE
     Write(outdat,*) '            Found ',M_FOUND,' Eigenvalues'
     Write(outdat,*)
     Call Print_Matrix(type_real_vector,W(1:M_Found),frmt='fr')
     write(vector_out) W(1:M_FOUND)           
  END IF
  IF (Get_Eigenvectors == .true.) THEN
      Call cpu_time(t_i)
      call DSTEIN(M_Size, D, E, Number_of_Eigenvectors , W, IBLOCK, ISPLIT, Z, M_SIze, &
                                WORK, IWORK, IFAIL, INFO)
      title = '                   Eigenvectors:'
      Call Print_Matrix(type_real_matrix,Z,M_Size, Number_of_Eigenvectors)
      write(vector_out) Z(:,1:Number_of_Eigenvectors)     
      Call cpu_time(t_f)
      write(outdat,*) '           ELapsed Time for eigenvectors= ',t_f-t_i      
  END IF
  DEALLOCATE( IBLOCK, ISPLIT, WORK, IWORK, IFAIL ) 
!********************************************************************************
  END SUBROUTINE H_0
!********************************************************************************
!********************************************************************************
  SUBROUTINE Propagation_CN_Length (cn_len)
  IMPLICIT NONE
  TYPE(CN_LENGTH)                            :: cn_len
  INTEGER                                    :: i
  INTEGER                                    :: INFO
  ALLOCATE (cn_len%ham_l(1:M_Size), cn_len%ham_u(1:M_Size), cn_len%ham_d(1:M_Size), cn_len%sol(1:M_Size ))

  pulse_time  =  delta_t * M_Size * 0.5d0
 
  cn_len%sol(1:m_size)    = Z (1:m_size,1)    
  cn_len%ham_l(1:m_size)  = E (1:m_size)
  cn_len%ham_u(1: m_size) = E (1:m_size)
  
  call cpu_time (t_i)
  time = 0.0d0
  DO i = 1, NO_Time_Steps
     IF (time < Pulse_Time) THEN
         E_t = 0.0d0
         IF (pulse == 'SQUARE' .OR. pulse == 'Square' .OR. pulse == 'square') THEN
             E_t   = E_0
         ELSE IF (pulse == 'SMOOTH' .OR. pulse == 'Smooth' .OR. pulse == 'smooth') THEN
             E_t   = E_0 * sin(PI * time / Pulse_Time)
         END IF
     END IF   
     cn_len%ham_d(1:m_size) = D(1:m_size) - x(1:m_size) * E_t *sin(Omega*(time + delta_t*0.5d0) + phase)
     time  = time + delta_t
     CALL CN(cn_len%ham_l,cn_len%ham_u,cn_len%ham_d,cn_len%sol,INFO)
  END DO

  call cpu_time (t_f)
  write (outcn,*) '***** Length Solution *****'
  IF ( t_f - t_i < 60 ) THEN
     write (outcn,*) 'Time elapsed: ', (t_f - t_i),'seconds'
  ELSE
     write (outcn,*) 'Time elapsed: ', (t_f - t_i)/real(60) , 'minutes'
  END IF
 
  IF (INFO .EQ. 0) THEN
  
!------------------------------ TEST 1 -----------------------------------!
!------------------ COMPUTING < LENGTH | H | LENGTH > --------------------!
!----IF AT E = 0, ENERGY IS EQUAL TO THE GROUND STATE ENERGY, THEN THE----!
!---------------------------TEST IS SUCCESSFULL---------------------------! 
  call Test_Energy_Real (cn_len%sol)
  write (outcn,*)
  write (outcn,*) '"TEST 1"' 
  write (outcn,*) 'Energy AT E = ', E_0, ': ', Energy
!-------------------------------------------------------------------------!
!-------------------------------------------------------------------------!
  title   = 'CN Length Gauge Solution Vector:'
  call prntcm(title, cn_len%sol, M_Size, 1, M_Size, 1, outcn)
  
  ELSE 
      print *, 'INFO is not zero in CN Length Gauge Solution'
  END IF   

  DEALLOCATE (cn_len%ham_l, cn_len%ham_u, cn_len%ham_d,cn_len%sol)
!********************************************************************************
  END SUBROUTINE Propagation_CN_Length 
!********************************************************************************
  SUBROUTINE Propagation_CN_Velocity (cn_vel)
  IMPLICIT NONE
  TYPE(CN_VELOCITY)                             :: cn_vel
  REAL(idp)                                     :: A
  REAL(idp), DIMENSION(:),         ALLOCATABLE  :: pt
  REAL(idp), DIMENSION(:),         ALLOCATABLE  :: wt
  REAL(idp), DIMENSION(:),         ALLOCATABLE  :: WORK
  COMPLEX(idp),DIMENSION(:)                     :: Velocity
  REAL(idp)                                     :: A_old
  INTEGER                                       :: i
  INTEGER                                       :: INFO

  Pulse_Time              =  delta_t * M_Size * 0.5d0
  Quad_Size = 4
  ALLOCATE ( HAM_D         (1 : M_Size),           &
             HAM_L         (1 : M_Size),           &
             HAM_U         (1 : M_Size),           &
             pt            (1 : Quad_Size),        &
             wt            (1 : Quad_Size),        &
             WORK          (1 : Quad_Size) )

  call cpu_time (t_0)
  
!  Velocity (: , 1)      = Z(: , 1)
  A = 0.0d0
  CALL Comp_A(A,0.0d0, delta_t*0.5d0, pt, wt,WORK)
  Velocity(:)           = Z(:,1)
  time = 0.0d0
  DO i =1, No_time_steps
     A_old = A
     HAM_D   (1 : M_Size)    = D(1 : M_Size) + A*A*0.5d0
     HAM_L   (1 : M_Size)    = E(1 : M_Size) - eye * A * 0.5d0 / Step_Size
     HAM_U   (1 : M_Size)    = E(1 : M_Size) + eye * A * 0.5d0 / Step_Size
     CALL CN(HAM_L, HAM_U, HAM_D, Velocity,INFO)
     CALL Comp_A (A, time + delta_t*0.5d0, time + delta_t*1.5d0, pt, wt, WORK)
     time  = time + delta_t
  END DO
  
  call cpu_time (t_f)
  write (outcn,*) '***** Velocity Solution *****'
  IF ( t_f - t_0 < 60 ) THEN
     write (outcn,*) 'Time elapsed: ', (t_f - t_0),'seconds'
  ELSE
     write (outcn,*) 'Time elapsed: ', (t_f - t_0)/real(60),'minutes'
  END IF

  IF (INFO .EQ. 0) THEN
!------------------------------ TEST 1 -----------------------------------!
!---------------- COMPUTING < VELOCITY | H | VELOCITY > ------------------!
!----IF AT E = 0, ENERGY IS EQUAL TO THE GROUND STATE ENERGY, THEN THE----!
!---------------------------TEST IS SUCCESSFULL---------------------------!

  call Test_Energy_Real (Velocity)
  write (outcn,*)
  write (outcn,*) '"TEST 1"'
  write (outcn,*) 'Energy AT E = ', E_0, ': ', Energy

!-------------------------------------------------------------------------!

  title   = ' Solution Vector (Velocity):'
  call prntcm(title, Velocity, M_Size, 1, M_Size, 1, outcn)

  ELSE 
      write(outcn,*) 'INFO is not zero in Velocity Gauge Solution'
  END IF 
  
  DEALLOCATE (WORK, wt, pt, HAM_D, HAM_U, HAM_L)
  END SUBROUTINE  Propagation_CN_Velocity
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
