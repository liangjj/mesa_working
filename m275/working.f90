!
  MODULE TDSE1D
  USE accuracy
  USE Matrix_Print
  IMPLICIT NONE
  COMPLEX(idp), PARAMETER                  :: eye = (0.0d0,1.0d0)
  REAL(idp),    PARAMETER                  :: PI = 3.141592653589793238462643383276D0
  INTEGER                                  :: M_Size       ! SIZE OF MATRIX
  REAL(idp)                                :: Step_Size
  LOGICAL                                  :: Get_Eigenvectors       ! GET EIGENVECTORS (TRUE)  
  LOGICAL                                  :: Print_Eigenvectors
  LOGICAL                                  :: Print_Solution
  INTEGER                                  :: indat=8         
  INTEGER                                  :: outdat=9          
  INTEGER                                  :: outcn = 10
  INTEGER                                  :: outso = 11     
  INTEGER                                  :: Number_of_Eigenvectors 
  REAL(idp), DIMENSION(:), ALLOCATABLE     :: X_N     
  REAL(idp), DIMENSION(:), ALLOCATABLE     :: D 
  REAL(idp), DIMENSION(:), ALLOCATABLE     :: E 
  REAL(idp), DIMENSION(:), ALLOCATABLE     :: W
  REAL(idp), DIMENSION(:), ALLOCATABLE     :: V
  REAL(idp), DIMENSION(:,:), ALLOCATABLE   :: Z 
  CHARACTER (LEN = 1)                      :: Y_N = 'N'
  INTEGER                                  :: vector_out
  INTEGER                                  :: IL
  INTEGER                                  :: IU
  INTEGER                                  :: M_Eig
  COMPLEX(idp),DIMENSION(:), ALLOCATABLE   :: HAM_U
  COMPLEX(idp),DIMENSION(:), ALLOCATABLE   :: HAM_L
  REAL(idp),   DIMENSION(:), ALLOCATABLE   :: HAM_D
  REAL(idp),   DIMENSION(:), ALLOCATABLE   :: Vector 
  REAL(idp)                                :: E_0
  REAL(idp)                                :: E_T 
  INTEGER                                  :: NO_Time_Steps
  REAL(idp)                                :: time
  CHARACTER(LEN=80)                        :: Method
  CHARACTER(LEN=80)                        :: pulse
  CHARACTER (LEN=80)                       :: title
  REAL(idp)                                :: Energy
  INTEGER                                  :: Quad_Size
  INTEGER                                  :: st
  REAL(idp)                                :: delta_t
  REAL(idp)                                :: Omega
  REAL(idp)                                :: phase
  REAL(idp)                                :: Pulse_Time         ! = DELTA_T * SIZE/2
  REAL(idp)                                :: t_0                ! INITIAL TIME (TIMER)
  REAL(idp)                                :: t_f                ! FINAL TIME (TIMER) 
  REAL(idp)                                :: left_end
  REAL(idp)                                :: right_end
  

!
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
  REAL(idp)                                :: VL
  REAL(idp)                                :: VU
  REAL(idp)                                :: DLAMCH
  REAL(idp)                                :: XX
  REAL(idp)                                :: AD
  REAL(idp)                                :: ABSTOL
  REAL                                     :: t_i
  REAL                                     :: t_f
  INTEGER                                  :: i 
  INTEGER                                  :: M_FOUND
  INTEGER                                  :: NSPLIT
  INTEGER                                  :: INFO
  CHARACTER (LEN = 1)                      :: Y_N = 'N'
  INTEGER                                  :: N
  INTEGER, DIMENSION(:), ALLOCATABLE       :: IBLOCK
  INTEGER, DIMENSION(:), ALLOCATABLE       :: ISPLIT
  TYPE(Real_Vector)                        :: type_real_vector
  TYPE(Real_Matrix)                        :: type_real_matrix
!
  open (vector_out,file='EigenVectors',access='sequential',form='unformatted',iostat=st,status='unknown')
  M_Size = (right_end - left_end)/step_size - 1
  ALLOCATE(D(1:M_Size), E(1:M_Size), X_N(1:M_Size), IBLOCK(1:M_Size), ISPLIT(1:M_Size),    &
           WORK(1:5*M_Size), W(1:M_Size), IWORK(1:5*M_Size), IFAIL(1:M_Size) )
  IF (Get_Eigenvectors == .true.) THEN
      Y_N = 'V'
      ALLOCATE(Z(M_Size,Number_of_Eigenvectors))
  END IF
  ABSTOL = 2 * DLAMCH('S')
!  XX = -N*Step_Size
  AD = - 1.d0 / ( 2.d0 * Step_Size * Step_Size )
  E(1:M_Size) = AD
  XX = 0.0d0
  DO i = 1, M_Size
     X_N(i) = XX
     D(i) =  - 2.d0 * AD !- 1.d0 / Sqrt (  ( 1.d0 + X_N(i)*X_N(i) ) )
     XX = XX + Step_Size
  END DO
  Call cpu_time(t_i)
  IF (IL < 1 ) Then
      IL = 1
  END IF 
  IF (IU > M_Size) THEN
      IU = M_Size
  END IF
  Call DSTEBZ('I', 'E', M_Size, VL, VU, IL, IU, ABSTOL, D, E, M_FOUND, NSPLIT, W,    &
               IBLOCK, ISPLIT, WORK, IWORK, INFO)
!  Call DSTEVX('N','A', M_Size, D, E, VL, VU, IL, IU, ABSTOL, M_FOUND, W, Z, M_Size,    &
!               WORK, IWORK, IFAIL, INFO)
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
      call DSTEIN(M_Size, D, E, Number_of_Eigenvectors , W, IBLOCK, ISPLIT, Z, M_Size, &
                                WORK, IWORK, IFAIL, INFO )
      IF (Print_Eigenvectors == .true. ) THEN
          title = '                   Eigenvectors:'
          Call Print_Matrix(type_real_matrix,Z,M_Size, Number_of_Eigenvectors)
      END IF
      write(vector_out) Z(:,1:Number_of_Eigenvectors)     
      Call cpu_time(t_f)
      write(outdat,*) '           ELapsed Time for eigenvectors= ',t_f-t_i      
  END IF
  DEALLOCATE( IBLOCK, ISPLIT, WORK, IWORK, IFAIL ) 
!********************************************************************************
  END SUBROUTINE H_0
!********************************************************************************
  SUBROUTINE Time_Dependent
  IMPLICIT NONE
  COMPLEX(idp),DIMENSION(:),   ALLOCATABLE :: Length            ! CN LENGTH SOLUTION
  REAL(idp),   DIMENSION(:),   ALLOCATABLE :: So_Real
  COMPLEX(idp),DIMENSION(:),   ALLOCATABLE :: S_Length
  COMPLEX(idp),DIMENSION(:),   ALLOCATABLE :: Velocity          ! CN VELOCITY SOLUTION
  COMPLEX(idp),DIMENSION(:),   ALLOCATABLE :: S_Velocity

  IF (Method == 'length') THEN
     Call H_0
     open(outcn,file='CN',access='sequential',form='formatted',iostat=st,status='unknown')
     ALLOCATE (Length(1 : M_Size))
     CALL Length_Solution (Length)
     DEALLOCATE (Length)
  ELSE IF (Method == 'velocity') THEN
     Call H_0
     open(outcn,file='CN',access='sequential',form='formatted',iostat=st,status='unknown')
     ALLOCATE (Velocity (1 : M_Size))
     CALL Velocity_Solution(Velocity)
     DEALLOCATE (Velocity)
  ELSE IF (Method == 'length_and_velocity') THEN
     Call H_0
     open(outcn,file='CN',access='sequential',form='formatted',iostat=st,status='unknown')
     ALLOCATE (Length(1 : M_Size), Velocity(1 : M_Size))
     CALL Length_Solution(Length)
     CALL Velocity_Solution(Velocity)
     DEALLOCATE (Length ,Velocity)
  ELSE IF (Method == 'real_time_split_operator') THEN
     CALL H_0
     open(outso,file='SO',access='sequential',form='formatted',iostat=st,status='unknown')
     ALLOCATE(S_Length(1 : M_Size))
     CALL REAL_Time_Split_Operator(S_Length)
     DEALLOCATE(S_Length)
  ELSE IF (Method == 'imag_time_split_operator') THEN
     CALL H_0
     open(outso,file='SO',access='sequential',form='formatted',iostat=st,status='unknown')
     ALLOCATE(So_Real(1 : M_Size))
     CALL IMAG_Time_Split_Operator(So_Real)
     DEALLOCATE(So_Real)
  ELSE IF (Method == 'diagonalize') THEN
       Call H_0
  ELSE
     CALL lnkerr('Method not available!')
  END IF
!********************************************************************************
  END SUBROUTINE Time_Dependent
!******************************************************************************** 
! ___________________________________________________________________________!
!_______________________________CRANK-NICOLSON_______________________________!
!                                                                            !
! If (first argument) FACT='N',the factored form of matrix A is not supplied on
! entry.
! If FACT='F',the factored form of matrix A is supplied on entry.
! If A*X=B (no transpose), then the second argument TRANS='N' 
! If A**T*X=B (transpose), then TRANS='T'
! If A**H*X=B (conjugate transpose), then TRANS='C'

  SUBROUTINE CN(HAM_L, HAM_U, HAM_D, X, INFO)
  IMPLICIT NONE
  REAL(idp),   DIMENSION(:)                      :: HAM_D
  COMPLEX(idp),DIMENSION(:)                      :: HAM_L
  COMPLEX(idp),DIMENSION(:)                      :: HAM_U
  COMPLEX(idp),DIMENSION(:), ALLOCATABLE         :: RHS_D    
  COMPLEX(idp),DIMENSION(:), ALLOCATABLE         :: RHS_L
  COMPLEX(idp),DIMENSION(:), ALLOCATABLE         :: RHS_U  

! ----------------------- ZGTSVX ARGUEMENTS ---------------------------- !
  CHARACTER                                      :: FACT
  CHARACTER                                      :: TRANS
  COMPLEX(idp),DIMENSION(:),   ALLOCATABLE       :: B
  COMPLEX(idp),DIMENSION(:),   ALLOCATABLE       :: D
  COMPLEX(idp),DIMENSION(:),   ALLOCATABLE       :: DF
  COMPLEX(idp),DIMENSION(:),   ALLOCATABLE       :: DL
  COMPLEX(idp),DIMENSION(:),   ALLOCATABLE       :: DLF
  COMPLEX(idp),DIMENSION(:),   ALLOCATABLE       :: DU
  COMPLEX(idp),DIMENSION(:),   ALLOCATABLE       :: DU2
  COMPLEX(idp),DIMENSION(:),   ALLOCATABLE       :: DUF
  COMPLEX(idp),DIMENSION(:),   ALLOCATABLE       :: WORK
  COMPLEX(idp),DIMENSION(:)                      :: X
  INTEGER                                        :: INFO
  INTEGER,     DIMENSION(:),   ALLOCATABLE       :: IPIV
  INTEGER                                        :: j
  INTEGER                                        :: NRHS
  REAL(idp),   DIMENSION(:),   ALLOCATABLE       :: BERR
  REAL(idp),   DIMENSION(:),   ALLOCATABLE       :: FERR
  REAL(idp)                                      :: RCOND
  REAL(idp),   DIMENSION(:),   ALLOCATABLE       :: RWORK

  NRHS = 1
  ALLOCATE(  RHS_L         (1 : M_Size),           &
             RHS_U         (1 : M_Size),           &    
             RHS_D         (1 : M_Size),           &
             D             (1 : M_Size),           &
             DL            (1 : M_Size),           &
             DU            (1 : M_Size),           &
             DLF           (1 : M_Size-1),         &
             DF            (1 : M_Size),           &
             DUF           (1 : M_Size),           &
             DU2           (1 : M_Size-2),         &
             IPIV          (1 : M_Size),           &
             B             (M_Size),               &
             FERR          (1 : NRHS),             &
             BERR          (1 : NRHS),             &
             WORK          (1 : 2*M_Size),         &
             RWORK         (1 : M_Size) )

  RHS_L   (1 : M_Size)    = -HAM_L(1 : M_Size) * eye * delta_t * 0.5d0
  RHS_U   (1 : M_Size)    = -HAM_U(1 : M_Size) * eye * delta_t * 0.5d0
  RHS_D   (1 : M_Size)    =  1 - HAM_D(1 : M_Size) * eye * delta_t * 0.5d0
  D       (1 : M_Size)    =  1 + HAM_D(1 : M_Size) * eye * delta_t * 0.5d0
  DL      (1 : M_Size)    =  HAM_L(1 : M_Size) * eye * delta_t * 0.5d0      
  DU      (1 : M_Size)    =  HAM_U(1 : M_Size) * eye * delta_t * 0.5d0

! .....................MULTIPLYING MATRIX AND VECTOR........................!     
  B(1)              = RHS_D(1)*X(1) + RHS_U(1)*X(2)
  DO j =2, M_Size-1
     B(j)           = RHS_L(j-1)*X(j-1) + RHS_D(j)*X(j) + RHS_U(j)*X(j+1)
  END DO
  B(M_Size)         = RHS_L(M_Size-1)*X(M_Size-1) + RHS_D(M_Size)*X(M_Size)    
!..........................................................................!
  CALL ZGTSVX('N','N', M_Size, NRHS, DL, D, DU, DLF, DF, DUF, DU2,     &
              IPIV, B, M_Size, X, M_Size, RCOND, FERR, BERR, WORK, RWORK, INFO)
    
  DEALLOCATE (RHS_L, RHS_U, RHS_D, D, DL, DU, DLF, DF, DUF, DU2, IPIV,  &
              B, FERR, BERR, WORK, RWORK)
!*******************************************************************************            
  END SUBROUTINE
!*******************************************************************************

!__________________________________________________________________________!
!__________________________ 2X2_MATRIX_EIGVAL_AND_EIGVEC___________________!
!__________________________________________________________________________!
!    GENERAL FORM: [a   b]                                                 !
!                  [c   d]                                                 !
!                                                                          !

  SUBROUTINE eigV(a,b,c,d, eigVal, eigVec)
  IMPLICIT NONE
  REAL(idp)                                      :: a
  REAL(idp)                                      :: b
  REAL(idp)                                      :: c
  REAL(idp)                                      :: d
  REAL(idp), DIMENSION(2)                        :: eigVal
  REAL(idp), DIMENSION(2,2)                      :: eigVec
  TYPE(Real_Vector)                              :: type_real_vector
  TYPE(Real_Matrix)                              :: type_real_matrix
  
  eigVec(1,1) = a
  eigVec(1,2) = b
  eigVec(2,1) = c
  eigVec(2,2) = d
!  Call Print_Matrix(type_real_matrix,eigVec,2,2,title='initial block')  
!  write (outdat,*) '****************************'
  eigVal(1)     = ( a+d + sqrt( (a+d)*(a+d) -4*(a*d - b*c) ) )* 0.5d0
  eigVal(2)     = ( a+d - sqrt( (a+d)*(a+d) -4*(a*d - b*c) ) )* 0.5d0

  eigVec(1,:)   = 1.0d0
  eigVec(2,:)   = ( eigVal(:) - a ) / b
  eigVec(:,1)   = 1.0d0/sqrt(eigVec(1,1)*eigVec(1,1) + eigVec(2,1)*eigVec(2,1))* eigVec(:,1)
  eigVec(:,2)   = 1.0d0/sqrt(eigVec(1,2)*eigVec(1,2) + eigVec(2,2)*eigVec(2,2))* eigVec(:,2)

!  write (outso,*)
!  Call Print_Matrix(type_real_vector,eigVal, frmt='fr', title='Eigevalues')  
!  write (outdat,*) '****************************'
!  Call Print_Matrix(type_real_matrix,eigVec,2,2,title='Eigenvectors')  
!  write (outdat,*) '****************************' 

!*******************************************************************************
  END SUBROUTINE eigV
!*******************************************************************************

!_____________________________2x2_BLOCK_MATRIX_____________________________!
!__________________________________________________________________________!
!                                                                          !
!    GENERAL FORM: [a   b]                                                 !
!                  [c   d]                                                 !
!                                                                          !

  SUBROUTINE REAL_TIME_BLOCK(Mtrx)
  IMPLICIT NONE
  INTEGER                                        :: n
  INTEGER                                        :: i
  INTEGER                                        :: j
  INTEGER                                        :: k
  REAL(idp), DIMENSION(:), ALLOCATABLE           :: eigVal
  REAL(idp), DIMENSION(:,:), ALLOCATABLE         :: eigVec
  COMPLEX(idp), DIMENSION(2,2, M_Size-1)         :: Mtrx
  CHARACTER (LEN=5)                              :: itoc 
 
  ALLOCATE (eigVal(2), eigVec(2,2))

  Mtrx(:,:,:)        = 0.0d0

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
  END SUBROUTINE REAL_TIME_BLOCK
!*******************************************************************************
  SUBROUTINE IMAG_TIME_BLOCK(Mtrx)
  IMPLICIT NONE
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
  END SUBROUTINE IMAG_TIME_BLOCK
!*******************************************************************************

!__________________________________________________________________________!
!________________________________TEST ENERGY_______________________________!
!_________________________COMPUTING  < X | H | X >_________________________!
!__________________________________________________________________________!
!                                                                          !
  
  SUBROUTINE Test_Energy_Real (Solution)
  IMPLICIT NONE
  COMPLEX(idp)                                    :: cdotc
  COMPLEX(idp),DIMENSION(:), ALLOCATABLE          :: MUL
  COMPLEX(idp),DIMENSION(:)                       :: Solution
  INTEGER                                         :: i
  
  ALLOCATE(MUL(M_Size))
  
  MUL(1)         = D(1) * Solution(1) + E(1) * Solution(2)
  Do i = 2, M_Size-1
     MUL(i)  = E(i-1) * Solution(i-1) + D(i) * Solution(i) + E(i+1) * Solution(i+1)
  END DO
  i = M_Size
  MUL(i)  = E(i-1) * Solution(i-1) + D(i ) * Solution(i) 
  Energy = cdotc(M_Size, Solution, 1, MUL, 1)
  DEALLOCATE(MUL)
!*******************************************************************************
  END SUBROUTINE Test_Energy_Real
!*******************************************************************************                                            
  SUBROUTINE Test_Energy_Complex (Solution)
  IMPLICIT NONE
  REAL(idp),DIMENSION(:), ALLOCATABLE          :: MUL
  REAL(idp),DIMENSION(:)                       :: Solution
  REAL(idp)                                    :: sdot
  INTEGER                                      :: i
  
  ALLOCATE(MUL(M_Size))
  
  MUL(1)         = D(1) * Solution(1) + E(1) * Solution(2)
  Do i = 2, M_Size-1
     MUL(i)  = E(i-1) * Solution(i-1) + D(i) * Solution(i) + E(i+1) * Solution(i+1)
  END DO
  i = M_Size
  MUL(i)  = E(i-1) * Solution(i-1) + D(i ) * Solution(i) 
  Energy = sdot(M_Size, Solution, 1, MUL, 1)
  DEALLOCATE(MUL)
!*******************************************************************************      
  END SUBROUTINE Test_Energy_Complex
!*******************************************************************************                                                                             
  SUBROUTINE Comp_A (A, t_l, t_u, pt, wt, WORK)
  IMPLICIT NONE
  INTEGER                                       :: i
  REAL(idp)                                     :: A
  REAL(idp)                                     :: alpha
  REAL(idp)                                     :: beta
  REAL(idp), DIMENSION(2)                       :: endpoints
  REAL(idp)                                     :: t_l
  REAL(idp)                                     :: t_u
  REAL(idp), DIMENSION(Quad_Size)               :: pt
  REAL(idp), DIMENSION(Quad_Size)               :: wt
  REAL(idp), DIMENSION(Quad_Size)               :: WORK

  CALL gaussq('legendre', Quad_Size, alpha, beta, 0, endpoints, WORK, pt, wt)
  CALL cnvtpt(pt, wt, t_l, t_u, Quad_size)
  A = 0.0d0
  DO i =1, Quad_Size
     E_t  =  0.0d0
     IF (pt(i) < Pulse_Time) THEN
            IF (pulse == 'SQUARE' .OR. pulse == 'Square' .OR. pulse == 'square') THEN
                E_t   = E_0
            ELSE IF (pulse == 'SMOOTH' .OR. pulse == 'Smooth' .OR. pulse =='smooth') THEN
                E_t   = E_0 * sin(PI * pt(i) / Pulse_Time)**2
            END IF
     END IF
     A = A - wt(i) *E_t*sin(Omega*pt(i) + phase)
  END DO
  
!********************************************************************************
  END SUBROUTINE Comp_A
!********************************************************************************
  SUBROUTINE Probability(PSI)
  IMPLICIT NONE
  COMPLEX(idp),DIMENSION(M_Size)                  :: PSI
  COMPLEX(idp)                                    :: cdotc
  INTEGER                                         :: i  
   
  DO i = 1, Number_of_Eigenvectors
     write(iout,1) i, abs(cdotc(M_Size,PSI,1,Z(1,i),1))**2
  END DO

  1 format('The probablity coefficient on state = ', i4, 2x, e15.8)
!********************************************************************************
  END SUBROUTINE
!********************************************************************************
  SUBROUTINE cnvtpt(pt,wt,t_l,t_u,Quad_size)
  IMPLICIT NONE
  INTEGER                                :: Quad_Size
  REAL(idp), DIMENSION(Quad_Size)        :: pt
  REAL(idp), DIMENSION(Quad_Size)        :: wt
  REAL(idp)                              :: t_l
  REAL(idp)                              :: t_u
  REAL(idp)                              :: f1
  REAL(idp)                              :: f2
  f1 = ( t_u - t_l )*.5D0
  f2 = ( t_l + t_u )*.5D0
  pt =  f1*pt + f2
  wt = wt*f1
!******************************************************************************** 
  END SUBROUTINE cnvtpt
!********************************************************************************
!__________________________________________________________________________!
!_______________________________LENGTH_SOLUTION____________________________!
!__________________________________________________________________________!
!                                                                          !

  SUBROUTINE Length_Solution (Length)
  IMPLICIT NONE
  INTEGER                                    :: i
  INTEGER                                    :: INFO
  COMPLEX(idp),DIMENSION(:)                  :: Length          ! CN Length Solution Vector (M_SIZE, NRHS)
  COMPLEX(idp),DIMENSION(:), ALLOCATABLE     :: New_Length      
              

  ALLOCATE (HAM_L(1 : M_Size), HAM_U(1 : M_Size), HAM_D(1 : M_Size), New_Length  (1 : M_Size ))

  Pulse_Time              =  delta_t * M_Size * 0.5d0
 
  Length  (:)             = Z (: , 1)    
  HAM_L   (1 : M_Size)    = E (1 : M_Size)
  HAM_U   (1 : M_Size)    = E (1 : M_Size)
  
  call cpu_time (t_0)
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
     HAM_D(:) = D(:) - x_n(:) * E_t *sin(Omega*(time + delta_t*0.5d0) + phase)
     time  = time + delta_t
     CALL CN(HAM_L, HAM_U, HAM_D, Length,INFO)
  END DO

  call cpu_time (t_f)
  write (outcn,*) '***** Length Solution *****'
  IF ( t_f - t_0 < 60 ) THEN
     write (outcn,*) 'Time elapsed: ', (t_f - t_0),'seconds'
  ELSE
     write (outcn,*) 'Time elapsed: ', (t_f - t_0)/real(60) , 'minutes'
  END IF
 
  IF (INFO .EQ. 0) THEN
  
!------------------------------ TEST 1 -----------------------------------!
!------------------ COMPUTING < LENGTH | H | LENGTH > --------------------!
!----IF AT E = 0, ENERGY IS EQUAL TO THE GROUND STATE ENERGY, THEN THE----!
!---------------------------TEST IS SUCCESSFULL---------------------------! 
  call Test_Energy_Real (Length)
  write (outcn,*)
  write (outcn,*) '"TEST 1"' 
  write (outcn,*) 'Energy AT E = ', E_0, ': ', Energy
!-------------------------------------------------------------------------!
!-------------------------------------------------------------------------!
  title   = 'CN Length Gauge Solution Vector:'
  call prntcm(title, Length, M_Size, 1, M_Size, 1, outcn)
  
  ELSE 
      print *, 'INFO is not zero in CN Length Gauge Solution'
  END IF   

  DEALLOCATE (HAM_L, HAM_U, HAM_D, New_Length)
!********************************************************************************
  END SUBROUTINE Length_Solution
!********************************************************************************
!_______________________________VELOCITY_SOLUTION__________________________!
!__________________________________________________________________________!

  SUBROUTINE velocity_Solution (Velocity)
  IMPLICIT NONE
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
  END SUBROUTINE Velocity_Solution
!********************************************************************************

!_______________________________ SPLIT OPERATOR ____________________________!
!___________________________________________________________________________!
!                                                                           !
!  LENGTH GAUGE:             
!      First_Block = [diag     off_d] 
!                    [off_d   diag/2]
!
!     Middle_Block = [diag/2  off_d] 
!                    [off_d  diag/2]        
!
!       Last_Block = [diag/2  off_d]
!                    [off_d    diag]     
!
!        [ v   0    0   ....]  
!     V =[ 0   v    0   ....] 
!        [ 0   0    v   ....]
!        [ ..  ..   ..  ....]      
!     v = - x_n * E * sin(omega*t + phase)   


  SUBROUTINE Real_Time_Split_Operator(Sol)
  IMPLICIT NONE
  INTEGER                                     :: i
  INTEGER                                     :: j
  INTEGER                                     :: k
  INTEGER                                     :: time_step
  COMPLEX(idp), DIMENSION(:)                  :: Sol
  COMPLEX(idp), DIMENSION(:,:,:), ALLOCATABLE :: Blocks         ! 2x2 even and odd matrices
  COMPLEX(idp), DIMENSION(:),     ALLOCATABLE :: tmp
  REAL(idp),    DIMENSION(:),     ALLOCATABLE :: V              ! = -x_n * E * sin(omega*t + phase) 
 
  
  ALLOCATE (Blocks (2,2,M_Size-1), tmp(1 : M_Size) ,V(1 : M_Size))

! LENGTH GAUGE
  
  CALL REAL_TIME_BLOCK(blocks)
  Sol(:)  = Z(:,1)
  DO j = 1, No_Time_steps
    !--------------------------------------------------------------------------!
     V(:)         = -x_n(:) * E_0 * sin(Omega*(time + delta_t *0.5d0) + phase)
     Sol(:)       = exp (-eye * V(:) * delta_t * 0.5d0) * Sol(:)
    !--------------------------------------------------------------------------!
    !-------------------------------ODD BLOCKS---------------------------------!    
     tmp(:)     = 0.0d0
     DO i = 1, M_Size - 2, 2
        tmp(i)    =  Blocks(1,1,i)*Sol(i) + Blocks(1,2,i)*Sol(i+1)
        tmp(i+1)  =  Blocks(2,1,i)*Sol(i) + Blocks(2,2,i)*Sol(i+1)
     END DO
     tmp(M_Size)        = Sol(M_Size)
    !-------------------------------------------------------------------------!
    !--------------------------------EVEN BLOCKS------------------------------!
   
     DO i = 2, M_Size - 1, 2
         Sol(i)    = Blocks(1,1,i)*tmp(i) + Blocks(1,2,i)*tmp(i+1)
         Sol(i+1)  = Blocks(2,1,i)*tmp(i) + Blocks(2,2,i)*tmp(i+1)
     END DO
     Sol(1)     =  tmp(1)
    !-------------------------------------------------------------------------!
    !-------------------------------ODD BLOCKS--------------------------------!
     DO i = 1, M_Size - 1, 2
         tmp(i)    =  Blocks(1,1,i)*Sol(i) + Blocks(1,2,i)*Sol(i+1)
         tmp(i+1)  =  Blocks(2,1,i)*Sol(i) + Blocks(2,2,i)*Sol(i+1)
     END DO
     tmp(M_Size)        = Sol(M_Size)
    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------! 
     Sol(:)       = exp(-eye * V(:) * delta_t * 0.5d0) * tmp(:)
    !-------------------------------------------------------------------------!
     time         = time + delta_t

  END DO
  call cpu_time (t_f)

  write (outso,*) ' (Real Time) Length Solution '
  IF ( t_f - t_0 < 60 ) THEN
     write (outso,*) 'Time elapsed: ', (t_f - t_0),'seconds'
  ELSE
     write (outso,*) 'Time elapsed: ', (t_f - t_0)/real(60) , 'minutes'
  END IF

!------------------------------ TEST 1 -----------------------------------!
!------------------ COMPUTING < LENGTH | H | LENGTH > --------------------!
!----IF AT E = 0, ENERGY IS EQUAL TO THE GROUND STATE ENERGY, THEN THE----!
!---------------------------TEST IS SUCCESSFULL---------------------------! 
  call Test_Energy_Real (Sol)
  write (outso,*)
  write (outso,*) '"TEST 1"'
  write (outso,*) 'Energy AT E = ', E_0, ': ', Energy
!------------------------------------------------------------------------!
  title   = ' Split method:'
  call prntcm(title, Sol, M_Size, 1, M_Size, 1, iout)

  DEALLOCATE (Blocks, tmp, V)
!********************************************************************************
  END SUBROUTINE Real_Time_Split_Operator
!********************************************************************************

  SUBROUTINE Imag_Time_Split_Operator(Sol)
  IMPLICIT NONE
  INTEGER                                     :: i
  INTEGER                                     :: j
  INTEGER                                     :: k
  INTEGER                                     :: time_step
  REAL(idp), DIMENSION(:)                     :: Sol
  REAL(idp), DIMENSION(:,:,:), ALLOCATABLE    :: Blocks         ! 2x2 even and odd matrices
  REAL(idp), DIMENSION(:),     ALLOCATABLE    :: tmp
  REAL(idp),    DIMENSION(:),     ALLOCATABLE :: V              ! = -x_n * E * sin(omega*t + phase) 
  REAL(idp)                                   :: old_Energy   
  CHARACTER(LEN=5)                            :: itoc        

  ALLOCATE (Blocks (2,2,M_Size-1), tmp(1 : M_Size) ,V(1 : M_Size))

! LENGTH GAUGE
  
  CALL IMAG_TIME_BLOCK(blocks)
  
!  DO i = 1, M_Size-1
!     title = 'Block = ' //itoc(i) 
!     CALL prntcm(title, blocks(1,1,i),2,2,2,2,outso)
!  END DO

  Sol(:)          = Z(:,1)
  CALL Test_Energy_Complex(Sol)
  old_Energy = Energy
  write(outso,*) 'Starting Energy: ', old_energy
  time = 0.0d0
  CALL cpu_time(t_0)
  DO j = 1, No_time_steps

    !--------------------------------------------------------------------------!
     V(:)         = -x_n(:) * E_0 * sin(Omega*(time + delta_t *0.5d0) + phase)
     Sol(:)       = exp (- V(:) * delta_t * 0.5d0) * Sol(:)
    !--------------------------------------------------------------------------!
    !-------------------------------ODD BLOCKS---------------------------------!    
     tmp(:)     = 0.0d0
     DO i = 1, M_Size - 2, 2
        tmp(i)    =  Blocks(1,1,i)*Sol(i) + Blocks(1,2,i)*Sol(i+1)
        tmp(i+1)  =  Blocks(2,1,i)*Sol(i) + Blocks(2,2,i)*Sol(i+1)
     END DO
     tmp(M_Size)        =  Sol(M_Size)
    !-------------------------------------------------------------------------!
    !--------------------------------EVEN BLOCKS------------------------------!
   
     DO i = 2, M_Size - 1, 2
         Sol(i)    = Blocks(1,1,i)*tmp(i) + Blocks(1,2,i)*tmp(i+1)
         Sol(i+1)  = Blocks(2,1,i)*tmp(i) + Blocks(2,2,i)*tmp(i+1)
     END DO
     Sol(1)  = tmp(1)
    !-------------------------------------------------------------------------!
    !-------------------------------ODD BLOCKS--------------------------------!
     DO i = 1, M_Size - 1, 2
         tmp(i)    =  Blocks(1,1,i)*Sol(i) + Blocks(1,2,i)*Sol(i+1)
         tmp(i+1)  =  Blocks(2,1,i)*Sol(i) + Blocks(2,2,i)*Sol(i+1)
     END DO
     tmp(M_Size)        = Sol(M_Size)
    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------! 
     Sol(:)       = exp(- V(:) * delta_t * 0.5d0) * tmp(:)
    !-------------------------------------------------------------------------!
     time         = time + delta_t

     CALL Test_Energy_Complex(Sol)
     write (outso,*) 'Step= ', j
     write (outso,*) 'Energy= ', Energy
     old_Energy = abs(old_Energy - Energy)
     write(outso,*) 'Absolute Energy Difference ', old_Energy
     IF (old_Energy <= 1.d-10) EXIT    
     old_Energy = Energy
  END DO
  call cpu_time (t_f)


!  IF ( t_f - t_0 < 60 ) THEN
!     write (outso,*) 'Time elapsed: ', (t_f - t_0),'seconds'
!  ELSE
!     write (outso,*) 'Time elapsed: ', (t_f - t_0)/real(60) , 'minutes'
!  END IF

!------------------------------ TEST 1 -----------------------------------!
!------------------ COMPUTING < LENGTH | H | LENGTH > --------------------!
!----IF AT E = 0, ENERGY IS EQUAL TO THE GROUND STATE ENERGY, THEN THE----!
!---------------------------TEST IS SUCCESSFULL---------------------------! 
  call Test_Energy_Complex (Sol)
  write (outso,*)
  write (outso,*) '"TEST 1"'
  write (outso,*) 'Energy AT E = ', E_0, ': ', Energy
!------------------------------------------------------------------------!
 
  DEALLOCATE (Blocks, tmp, V)
!********************************************************************************
  END SUBROUTINE IMAG_Time_Split_Operator
!********************************************************************************
  END MODULE TDSE1D
!********************************************************************************
