!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! full_matrix_vector_multiply_module
!**begin prologue     full_matrix_vector_multiply_module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Contains the routines needed to multiply the
!***                  Hamiltonian on a vector or to mutiply the exponential
!***                  sector propagators on a vector. Explicit interfaces are 
!***                  used to allow a transparent use of generic subroutines 
!***                  which work for both real and complex vectors.  
!***                  This feature permits a single code to be used for both 
!***                  real and imaginary time propagation.
!***description       See the specific routined.
!**references
!**modules needed     See USE statements below
!**end prologue       full_matrix_vector_multiply_module
!***********************************************************************
!***********************************************************************
                      MODULE full_matrix_vector_multiply_module
                         USE dvrprop_global
                         USE dvr_shared
                         USE dvr_global
                         USE Iterative_Global
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE matrix_vector_multiply
              MODULE PROCEDURE full_matrix_on_vector_d,                &
                               full_matrix_on_vector_z,                &
                               full_matrix_on_vector_zz,               &
                               lower_triangle_matrix_on_vector_d,      &
                               lower_triangle_matrix_on_vector_z,      &   
                               lower_triangle_matrix_on_vector_zz   
                END INTERFACE matrix_vector_multiply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck full_matrix_on_vector_d
!***begin prologue     full_matrix_on_vector_d    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine is the driver routine which 
!***                   calculates,
!------------------------------------------------------------------------------------
!
!                 V = H * V
!
!------------------------------------------------------------------------------------
!***                   Where H is a general matrix.
!***                   The parameter nj is a dummy and can take on values
!***                   consistent with any dimensional problem.  In a 1D
!***                   case nj=1, in 2D nj=nx and in 3D nj=ny*nx.
!***                   It is assumed that double counting of the diagonals
!***                   has been accounted for by halving the elements at
!***                   the interior, connecting, DVR element points.
!***references
!***routines called    
!
!***end prologue       full_matrix_on_vector_d                     
!
  SUBROUTINE full_matrix_on_vector_d (matrix,vector_in,vector_out)
  IMPLICIT NONE
  REAL*8, DIMENSION(:,:)                  :: matrix
  REAL*8, DIMENSION(:)                    :: vector_in
  REAL*8, DIMENSION(:)                    :: vector_out
!
!
  CALL ebc(vector_out,matrix,vector_in,n3d,n3d,1)
!
 END SUBROUTINE full_matrix_on_vector_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!deck full_matrix_on_vector_z
!***begin prologue     full_matrix_on_vector_z    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine is the driver routine which 
!***                   calculates,
!------------------------------------------------------------------------------------
!
!                 V = H * V
!
!------------------------------------------------------------------------------------
!***                   Where H is a general matrix.
!
!***                   The parameter nj is a dummy and can take on values
!***                   consistent with any dimensional problem.  In a 1D
!***                   case nj=1, in 2D nj=nx and in 3D nj=ny*nx.
!***                   It is assumed that double counting of the diagonals
!***                   has been accounted for by halving the elements at
!***                   the interior, connecting, DVR element points.
!***references
!***routines called    
!
!***end prologue       full_matrix_on_vector_z                     

  SUBROUTINE full_matrix_on_vector_z(matrix,vector_in,vector_out)
  IMPLICIT NONE
  REAL*8,     DIMENSION(:,:)            :: matrix  
  COMPLEX*16, DIMENSION(:)              :: vector_in
  COMPLEX*16, DIMENSION(:)              :: vector_out
!
!
  CALL ebcc(vector_out,matrix,vector_in,n3d,n3d,1) 
!
  END SUBROUTINE full_matrix_on_vector_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!deck full_matrix_on_vector_zz
!***begin prologue     full_matrix_on_vector_zz
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine is the driver routine which 
!***                   calculates,
!------------------------------------------------------------------------------------
!
!                 V = H * V
!
!------------------------------------------------------------------------------------
!***                   Where H is a general matrix.
!
!***                   The parameter nj is a dummy and can take on values
!***                   consistent with any dimensional problem.  In a 1D
!***                   case nj=1, in 2D nj=nx and in 3D nj=ny*nx.
!***                   It is assumed that double counting of the diagonals
!***                   has been accounted for by halving the elements at
!***                   the interior, connecting, DVR element points.
!***references
!***routines called    
!
!***end prologue       full_matrix_on_vector_zz
!
  SUBROUTINE full_matrix_on_vector_zz(matrix,vector_in,vector_out)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:,:)            :: matrix  
  COMPLEX*16, DIMENSION(:)              :: vector_in
  COMPLEX*16, DIMENSION(:)              :: vector_out
!
!
  CALL cebc(vector_out,matrix,vector_in,n3d,n3d,1) 
!
  END SUBROUTINE full_matrix_on_vector_zz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck lower_triangle_matrix_on_vector_d
!***begin prologue     lower_triangle_matrix_on_vector_d    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine is the driver routine which 
!***                   calculates,
!------------------------------------------------------------------------------------
!
!                 V = H * V
!
!------------------------------------------------------------------------------------
!***                   Where H is a real, symmetric matrix with only the triangle stored.
!***references
!***routines called    
!
!***end prologue       lower_triangle_matrix_on_vector_d                     
!
  SUBROUTINE lower_triangle_matrix_on_vector_d (lower_triangle_matrix,vector_in,     &
                                                vector_out)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                    :: lower_triangle_matrix
  REAL*8, DIMENSION(:)                    :: vector_in
  REAL*8, DIMENSION(:)                    :: vector_out
  INTEGER                                 :: i
  INTEGER                                 :: j
  INTEGER                                 :: count
!
!
  vector_out(:) = 0.d0
  count = 0
  DO i=1,n3d
     DO j=1,i-1
        count = count + 1
        vector_out(i) = vector_out(i) + lower_triangle_matrix(count) * vector_in(j)
        vector_out(j) = vector_out(j) + lower_triangle_matrix(count) * vector_in(i)
     END DO
     count = count + 1
     vector_out(i) = vector_out(i) + lower_triangle_matrix(count) * vector_in(i)     
  END DO
!
 END SUBROUTINE lower_triangle_matrix_on_vector_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck lower_triangle_matrix_on_vector_z
!***begin prologue     lower_triangle_matrix_on_vector_z    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine is the driver routine which 
!***                   calculates,
!------------------------------------------------------------------------------------
!
!                 V = H * V
!
!------------------------------------------------------------------------------------
!***                   Where H is a symmetric matrix with only the triangle stored.
!***references
!***routines called    
!
!***end prologue       lower_triangle_matrix_on_vector_z                     
!
  SUBROUTINE lower_triangle_matrix_on_vector_z (lower_triangle_matrix,vector_in,     &
                                                vector_out)
  IMPLICIT NONE
  REAL*8,     DIMENSION(:)                :: lower_triangle_matrix
  COMPLEX*16, DIMENSION(:)                :: vector_in
  COMPLEX*16, DIMENSION(:)                :: vector_out
  INTEGER                                 :: i
  INTEGER                                 :: j
  INTEGER                                 :: count
!
!
  vector_out(:) = (0.d0,0.d0)
  count = 0
  DO i=1,n3d
     DO j=1,i-1
        count = count + 1
        vector_out(i) = vector_out(i) + lower_triangle_matrix(count)        &
                                      *                                     &
                                        vector_in(j)
        vector_out(j) = vector_out(j) + lower_triangle_matrix(count)        &
                                      *                                     &
                                        vector_in(i)
     END DO
     count = count + 1
     vector_out(i) = vector_out(i) + lower_triangle_matrix(count) * vector_in(i)     
  END DO
!
 END SUBROUTINE lower_triangle_matrix_on_vector_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck lower_triangle_matrix_on_vector_zz
!***begin prologue     lower_triangle_matrix_on_vector_zz
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine is the driver routine which 
!***                   calculates,
!------------------------------------------------------------------------------------
!
!                 V = H * V
!
!------------------------------------------------------------------------------------
!***                   Where H is a hermitian matrix with only the triangle stored.
!***references
!***routines called    
!
!***end prologue       lower_triangle_matrix_on_vector_zz                     
!
  SUBROUTINE lower_triangle_matrix_on_vector_zz (lower_triangle_matrix,vector_in,     &
                                                 vector_out)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)                :: lower_triangle_matrix
  COMPLEX*16, DIMENSION(:)                :: vector_in
  COMPLEX*16, DIMENSION(:)                :: vector_out
  INTEGER                                 :: i
  INTEGER                                 :: j
  INTEGER                                 :: count
!
!
  vector_out(:) = (0.d0,0.d0)
  count = 0
  DO i=1,n3d
     DO j=1,i-1
        count = count + 1
        vector_out(i) = vector_out(i) + lower_triangle_matrix(count)        &
                                      *                                     &
                                        vector_in(j)
        vector_out(j) = vector_out(j) + conjg(lower_triangle_matrix(count)) &
                                      *                                     &
                                        vector_in(i)
     END DO
     count = count + 1
     vector_out(i) = vector_out(i) + lower_triangle_matrix(count) * vector_in(i)     
  END DO
!
 END SUBROUTINE lower_triangle_matrix_on_vector_zz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            END       MODULE full_matrix_vector_multiply_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
