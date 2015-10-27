!***********************************************************************
! Matrix_Manipulation_Module
!**begin prologue     Matrix_Manipulation_Module
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
!**end prologue       Matrix_Manipulation_Module
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   MODULE Matrix_Manipulation_Module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              INTERFACE vexp
              MODULE PROCEDURE vexp_d,                                      &
                               vexp_z
              END INTERFACE  vexp
              INTERFACE mvmul
              MODULE PROCEDURE mvmul_d,                                     &
                               mvmul_z,                                     &
                               mvmul_zz                                      
              END INTERFACE  mvmul
              INTERFACE ebct
              MODULE PROCEDURE ebct_d,                                      &
                               ebct_z
              END INTERFACE  ebct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ebct_d(a,b,c,ni,nk,nj)
!***begin prologue     ebct_d
!***date written       850601  (yymmdd)
!***revision date      yymmdd  (yymmdd)
!***keywords           matrix, multiply,
!***author             saxe, paul (lanl)
!***source             @(#)ebct.f	5.1   11/6/94
!***purpose                                               t
!                      vectorized matrix operation:  a=b*c .
!***description
!                      call ebct_d(a,b,c,ni,nk,nj)
!                        a       output matrix, (ni,nj).
!                        b       input matrix, (ni,nk).
!                        c       input matrix, (nj,nk).
!
!***references
!***routines called    dgemm(clams)
!***end prologue       ebct_d
  implicit none
!
  integer                       :: ni, nj, nk
  real*8, dimension(ni,nj)      :: a
  real*8, dimension(ni,nk)      :: b
  real*8, dimension(nj,nk)      :: c
  real*8                        :: zero=0.d0, one=1.d0
!
!
  call dgemm('n','t',ni,nj,nk,one,b,ni,c,nj,zero,a,ni)
!
END SUBROUTINE ebct_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck ebct_z
!***begin prologue     ebct_z
!***date written       880423   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           complex matrix multiply
!***author             schneider, barry (lanl)
!***source             mylib
!                                                t
!***purpose            matrix multiply  a = b * c
!***description        vectorized matrix multiply
!***                   for complex a, b and c
!***                  
!***references         none
!
!***routines called    zgemm
!***end prologue       ebct_z
  subroutine ebct_z(a,b,c,ni,nk,nj)
  implicit none
  integer                           :: ni, nj, nk
  complex*16, dimension(ni,nj)      :: a
  complex*16, dimension(ni,nk)      :: b
  complex*16, dimension(nj,nk)      :: c
  complex*16                        :: zero=(0.d0,0.d0), one=(1.d0,0.d0)
  call zgemm('n','t',ni,nj,nk,one,b,ni,c,nj,zero,a,ni)
END SUBROUTINE ebct_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  subroutine mvmul_d(v,matin,matout,n,m)
!***begin prologue     mvmul_d
!***date written       850601  (yymmdd)
!***revision date      yymmdd  (yymmdd)
!***keywords           vector, multiply
!***author             schneider, barry (nsf)
!***source             @(#)mvmul.f	1.1   6/22/93
!***purpose            vectorized matrix vector multiply:  matout=matin*v.
!***description
!                      call mvmul(v,matin,matout,n,m)
!                        v        input vector of length m.
!                        matin    input matrix of size (n*m).
!                        matout   output matrix of size (n*m).
!                        m        length of vector ( matrix column size ).
!                        n        matrix row size. 
!c
!***references
!***routines called    (none)
!***end prologue       mvmul_d
  REAL*8, DIMENSION(m)              :: v
  REAL*8, DIMENSION(n,m)            :: matin
  REAL*8, DIMENSION(n,m)            :: matout
  INTEGER                           :: i
  DO i=1,n
     matout(i,:) = matin(i,:) * v(:)
  END DO
END SUBROUTINE mvmul_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  subroutine mvmul_z(v,matin,matout,n,m)
!***begin prologue     mvmul_z
!***date written       850601  (yymmdd)
!***revision date      yymmdd  (yymmdd)
!***keywords           vector, multiply
!***author             schneider, barry (nsf)
!***source             @(#)mvmul.f	1.1   6/22/93
!***purpose            vectorized matrix vector multiply:  matout=matin*v.
!***description
!                      call mvmul(v,matin,matout,n,m)
!                        v        input vector of length m.
!                        matin    input matrix of size (n*m).
!                        matout   output matrix of size (n*m).
!                        m        length of vector ( matrix column size ).
!                        n        matrix row size. 
!c
!***references
!***routines called    (none)
!***end prologue       mvmul_z
  REAL*8,     DIMENSION(m)           :: v
  COMPLEX*16, DIMENSION(n,m)         :: matin
  COMPLEX*16, DIMENSION(n,m)         :: matout
  INTEGER                            :: i
  DO i=1,n
     matout(i,:) = matin(i,:) * v(:)
  END DO
END SUBROUTINE mvmul_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  subroutine mvmul_zz(v,matin,matout,n,m)
!***begin prologue     mvmul_zz
!***date written       850601  (yymmdd)
!***revision date      yymmdd  (yymmdd)
!***keywords           vector, multiply
!***author             schneider, barry (nsf)
!***source             @(#)mvmul.f	1.1   6/22/93
!***purpose            vectorized matrix vector multiply:  matout=matin*v.
!***description
!                      call mvmul(v,matin,matout,n,m)
!                        v        input vector of length m.
!                        matin    input matrix of size (n*m).
!                        matout   output matrix of size (n*m).
!                        m        length of vector ( matrix column size ).
!                        n        matrix row size. 
!c
!***references
!***routines called    (none)
!***end prologue       mvmul_zz
  COMPLEx*16, DIMENSION(m)           :: v
  COMPLEX*16, DIMENSION(n,m)         :: matin
  COMPLEX*16, DIMENSION(n,m)         :: matout
  INTEGER                            :: i
  DO i=1,n
     matout(i,:) = matin(i,:) * v(:)
  END DO
END SUBROUTINE mvmul_zz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   subroutine vexp_d(v,w,tau,n)
!***begin prologue     vexp_d
!***date written       850601  (yymmdd)
!***revision date      yymmdd  (yymmdd)
!***keywords           vector, exponentiation
!***author             saxe, paul (lanl)
!***source             @(#)vexp.f	5.1   11/6/94
!***purpose            vectorized vector exponentiation:  v=exp(w) .
!***description
!                      call vexp(v,w,tau,n)
!                        v        output vector of length (n).
!                        w        input vector of length (n).
!                        tau      scalar.
!                        n        vector length.
!
!***references
!***routines called    (none)
!***end prologue       vexp_d
  REAL*8, DIMENSION(n)                  :: v
  REAL*8, DIMENSION(n)                  :: w
  REAL*8                                :: tau
  v(:) = exp ( tau * w(:) )
  return
END SUBROUTINE vexp_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   subroutine vexp_z(v,w,tau,n)
!***begin prologue     vexp_z
!***date written       850601  (yymmdd)
!***revision date      yymmdd  (yymmdd)
!***keywords           vector, exponentiation
!***author             saxe, paul (lanl)
!***source             @(#)vexp.f	5.1   11/6/94
!***purpose            vectorized vector exponentiation:  v=exp(w) .
!***description
!                      call vexp(v,w,tau,n)
!                        v        output vector of length (n).
!                        w        input vector of length (n).
!                        tau      scalar.
!                        n        vector length.
!
!***references
!***routines called    (none)
!***end prologue       vexp_z
  COMPLEX*16, DIMENSION(n)              :: v
  REAL*8,     DIMENSION(n)              :: w
  REAL*8                                :: tau
  v(:) = exp ( tau * w(:) )
  return
END SUBROUTINE vexp_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END MODULE matrix_manipulation_module
