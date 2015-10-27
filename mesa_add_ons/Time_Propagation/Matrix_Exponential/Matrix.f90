  Program Matrix
  USE input_output
  USE Matrix_Exponential_Module
  IMPLICIT NONE
  INTEGER                                   :: n, i, j
  REAL*8, DIMENSION(:,:), ALLOCATABLE       :: A_Real, Exp_Real
  REAL*8                                    :: t
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE   :: A_Hermitian, Exp_Complex
  CHARACTER(LEN=16)                         :: type
  CHARACTER(LEN=80)                         :: title
  OPEN(inp,file='matrix.inp',status='old')
  OPEN(iout,file='matrix.out',status='unknown')
  READ(inp,*) type, n, t
  write(iout,*) 'Type Matrix = ',type,'Size Matrix = ',n, 'Time Step =',t
  IF (type == 'real') THEN
     ALLOCATE ( A_Real(n,n), Exp_Real(n,n) )
     DO i=1,n
        DO j=1,i
           A_Real(i,j) = 1.d0/(i+j)
           A_Real(j,i) = A_Real(i,j)
        END DO
     END DO   
     title='Real Matrix'
     call prntfmn(title,A_Real,n,n,n,n,iout,'e')
     CALL Matrix_Exponential(A_Real,Exp_Real,t,n)
     DEALLOCATE ( A_Real, Exp_Real )
  ELSE
     ALLOCATE ( A_Hermitian(n,n), Exp_Complex(n,n) )
     DO i=1,n
        DO j=1,i
           A_Hermitian(i,j) = 1.d0/(i+j)
           A_Hermitian(j,i) = A_Hermitian(i,j)
        END DO
     END DO
     title='Hermitian Matrix'
     call prntcmn(title,A_Hermitian,n,n,n,n,iout,'e')
     CALL Matrix_Exponential(A_Hermitian,Exp_Complex,t,n)
     DEALLOCATE ( A_Hermitian, Exp_Complex )
  END IF
  END Program Matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
