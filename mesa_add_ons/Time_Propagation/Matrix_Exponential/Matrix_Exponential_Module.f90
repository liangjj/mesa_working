!***********************************************************************
! Matrix_Exponential_Module
!**begin prologue     Matrix_Exponential_Module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            To form a matrix representation of the operator
!**
!**                         Exp( tau * A )
!**                   Note that tau may be equal to i * t, in which
!**                   case we are computing the complex exponential.
!***description       In the first step we diagonaliz A.  The we use the
!***                  the eigenvalues and eigenvectors to transform the
!***                  matrix to the original representation.
!**references
!**modules needed     See USE statements below
!**end prologue       regional_module
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   MODULE Matrix_Exponential_Module
                   USE Input_Output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              INTERFACE matrix_exponential
              MODULE PROCEDURE matrix_exponential_d,                        &
                               matrix_exponential_h
              END INTERFACE  matrix_exponential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        CONTAINS
!deck matrix_exponential_d
!***begin prologue     matrix_exponential_d    
!***date written       020506   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           exponential of a matrix
!***author             schneider, b, i. (nsf)
!***source             matrix_exponential
!***purpose            Form the exponential of a real symmetrix matrix
!**description         This routine computes the matrix representation of
!**                    
!**                         Exp ( - tau * A )
!**                    where A is a real symmetric matrix.  A is not destroyed by
!**                    this routine.
!***routines called    
!
!***end prologue       matrix_exponential_d
!
  SUBROUTINE matrix_exponential_d(A,Exp_A,tau,n)
  IMPLICIT NONE
  INTEGER                                :: n, i, info
  REAL*8, DIMENSION(:,:)                 :: A
  REAL*8, DIMENSION(:,:)                 :: Exp_A
  REAL*8, DIMENSION(:),   ALLOCATABLE    :: Eigen_Values
  REAL*8, DIMENSION(:,:), ALLOCATABLE    :: Eigen_Vectors
  REAL*8, DIMENSION(:,:), ALLOCATABLE    :: Temp
  REAL*8, DIMENSION(:),   ALLOCATABLE    :: Scratch
  REAL*8                                 :: tau
  CHARACTER (LEN=80)                     :: title
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         Step One is to digonalize A
!
! Allocate some storage for diagonalization routine.
!
  ALLOCATE ( Eigen_Values(n), Eigen_Vectors(n,n), Scratch(10*n), Temp(n,n) )  
  Eigen_Vectors = A
  CALL dsyev('v','l',n,Eigen_Vectors,n,Eigen_Values,           &
              Scratch,10*n,info)
!
! Print some information if you want to.  This can be deleted
!
  title='Eigenvalues'
  call prntfmn(title,Eigen_Values,n,1,n,1,iout,'e')
  title='Eigenvectors'
  call prntfmn(title,Eigen_Vectors,n,n,n,n,iout,'e')
!
! Form an intermediate quantity so that the final matrix multiplication
! can be easily done.
!
  DO i=1,n
     Temp(:,i) = exp ( -tau * Eigen_Values(i) ) * Eigen_Vectors(:,i)
  END DO
!
! Call the appropriate matrix multiplication routine
!
  CALL ebct(Exp_A,Eigen_Vectors,Temp,n,n,n)
!
! Deallocate the unneeded storage
!
  DEALLOCATE ( Eigen_Values, Eigen_Vectors, Scratch, Temp )  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE matrix_exponential_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck matrix_exponential_h
!***begin prologue     matrix_exponential_h    
!***date written       020506   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           exponential of a matrix
!***author             schneider, b, i. (nsf)
!***source             matrix_exponential
!***purpose            Form the exponential of a Hermitian matrix
!**description         This routine computes the matrix representation of
!**                    
!**                         Exp ( - i * tau * A )
!**                    where A is a Hermitianc matrix.  A is not destroyed by
!**                    this routine.
!***routines called    
!
!***end prologue       matrix_exponential_h
!
  SUBROUTINE matrix_exponential_h(A,Exp_A,t,n)
  IMPLICIT NONE
  INTEGER                                    :: n, i, info
  COMPLEX*16, DIMENSION(:,:)                 :: A
  COMPLEX*16, DIMENSION(:,:)                 :: Exp_A
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE    :: Eigen_Vectors
  REAL*8,     DIMENSION(:),   ALLOCATABLE    :: Eigen_Values
  COMPLEX*16, DIMENSION(:),   ALLOCATABLE    :: Work
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE    :: Temp
  REAL*8,     DIMENSION(:),   ALLOCATABLE    :: Rwork
  REAL*8                                     :: t
  CHARACTER (LEN=80)                         :: title
  COMPLEX*16                                 :: eye=(1.d0,0.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         Step One is to digonalize A
!
! Allocate some storage for diagonalization routine.
!
  ALLOCATE ( Eigen_Values(n), Eigen_Vectors(n,n), Work(10*n), Rwork(10*n), Temp(n,n) )  
  Eigen_Vectors = A
  CALL zheev('v','l',n,Eigen_Vectors,n,Eigen_Values,              &
              Work,10*n,Rwork,info)
!
! Print some information if you want to.  This can be deleted
!
  title='Eigenvalues'
  call prntfmn(title,Eigen_Values,n,1,n,1,iout,'e')
  title='Eigenvectors'
  call prntcmn(title,Eigen_Vectors,n,n,n,n,iout,'e')
!
! Form an intermediate quantity so that the final matrix multiplication
! can be easily done.
!
  DO i=1,n
     Temp(:,i) = exp ( - eye * t * Eigen_Values(i) ) * Eigen_Vectors(:,i)
  END DO
!
! Call the appropriate matrix multiplication routine
!
  CALL cebct(Exp_A,Eigen_Vectors,Temp,n,n,n)
!
! Deallocate the unneeded storage
!
  DEALLOCATE ( Eigen_Values, Eigen_Vectors, Work, Rwork, Temp )    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE matrix_exponential_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END MODULE matrix_exponential_module
