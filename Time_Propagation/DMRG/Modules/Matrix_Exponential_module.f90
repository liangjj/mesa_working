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
!***description       In the first step we diagonalize A.  Then we use the
!***                  the eigenvalues and eigenvectors to transform the
!***                  matrix to the original representation.
!**references
!**modules needed     See USE statements below
!**end prologue       regional_module
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   MODULE Matrix_Exponential_Module
                   USE Matrix_Manipulation_Module
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
  SUBROUTINE matrix_exponential_d(A,Exp_A,delta_t,n)
  IMPLICIT NONE
  INTEGER                                :: n, i, j, k, info
  REAL*8, DIMENSION(:,:)                 :: A
  REAL*8, DIMENSION(:,:)                 :: Exp_A
  REAL*8, DIMENSION(:),   ALLOCATABLE    :: Eigen_Values
  REAL*8, DIMENSION(:,:), ALLOCATABLE    :: Eigen_Vectors
  REAL*8, DIMENSION(:,:), ALLOCATABLE    :: Temp
  REAL*8, DIMENSION(:,:), ALLOCATABLE    :: Scratch
  REAL*8                                 :: delta_t
  REAL*8                                 :: tau
  INTEGER                                :: lwork
  CHARACTER (LEN=80)                     :: title
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         Step One is to digonalize A
!
! Allocate some storage for diagonalization routine.
!
  lwork = max(n*n,10*n)
  ALLOCATE ( Eigen_Values(n), Eigen_Vectors(n,n), Scratch(lwork,1) )  
  Eigen_Vectors = A
  CALL dsyev('v','l',n,Eigen_Vectors,n,Eigen_Values,Scratch,lwork,info)
!
! Print some information if you want to.  This can be deleted
!
!  title='Eigenvalues'
!  call prntfmn(title,Eigen_Values,n,1,n,1,iout,'e')
!  title='Eigenvectors'
!  call prntfmn(title,Eigen_Vectors,n,n,n,n,iout,'e')
!
! Form an intermediate quantity so that the final matrix multiplication
! can be easily done.
!
  tau = - delta_t
  CALL vexp(Eigen_Values,Eigen_Values,tau,n)
  CALL mvmul(Eigen_Values,Eigen_Vectors,Scratch,n,n)
  CALL ebct(Exp_A,Eigen_Vectors,Scratch,n,n,n)
!
  DEALLOCATE ( Eigen_Values, Eigen_Vectors, Scratch )  
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
  SUBROUTINE matrix_exponential_h(A,Exp_A,delta_t,n)
  IMPLICIT NONE
  INTEGER                                    :: n, i, j, k, info

  COMPLEX*16, DIMENSION(:,:)                 :: A
  COMPLEX*16, DIMENSION(:,:)                 :: Exp_A
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE    :: Eigen_Vectors
  REAL*8,     DIMENSION(:),   ALLOCATABLE    :: Eigen_Values
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE    :: Scratch
  COMPLEX*16, DIMENSION(:),   ALLOCATABLE    :: Zwork
  REAL*8,     DIMENSION(:),   ALLOCATABLE    :: Rwork
  REAL*8                                     :: delta_t
  REAL*8                                     :: tau
  CHARACTER (LEN=80)                         :: title
  COMPLEX*16                                 :: eye=(0.d0,1.d0)
  COMPLEX*16                                 :: expeig
  INTEGER                                    :: lwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         Step One is to digonalize A
!
! Allocate some storage for diagonalization routine.
!
  lwork=max(n*n,10*n)
  ALLOCATE ( Eigen_Values(n), Eigen_Vectors(n,n), Scratch(n,n), Zwork(n), Rwork(lwork) )  
  Eigen_Vectors = A
  CALL zheev('v','l',n,Eigen_Vectors,n,Eigen_Values,Scratch,lwork,Rwork,info)
!
  tau = - delta_t
  CALL vexp(Zwork,Eigen_Values,tau,n)
  CALL mvmul(Zwork,Eigen_Vectors,Scratch,n,n)
  CALL ebct(Exp_A,Eigen_Vectors,Scratch,n,n,n)
!
! Deallocate the unneeded storage
!
  DEALLOCATE ( Eigen_Values, Eigen_Vectors, Scratch, Zwork, Rwork )    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE matrix_exponential_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END MODULE matrix_exponential_module
