!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! matrix_vector_multiply_module
!**begin prologue     matrix_vector_multiply_module
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
!**end prologue       arnoldi_module
!***********************************************************************
!***********************************************************************
                      MODULE matrix_vector_multiply_module
                         USE dvrprop_global
                         USE dvr_shared
                         USE dvr_global
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE finite_element_m_v
              MODULE PROCEDURE finite_element_m_v_d,                  &
                               finite_element_m_v_z
                END INTERFACE finite_element_m_v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE dvr_mat_mul
              MODULE PROCEDURE dvr_mat_mul_d,                           &
                               dvr_mat_mul_z
                END INTERFACE dvr_mat_mul
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE v_mat_v_gen
              MODULE PROCEDURE v_mat_v_gen_d,                           &
                               v_mat_v_gen_z
                END INTERFACE  v_mat_v_gen
                    INTERFACE v_mat_v_2
              MODULE PROCEDURE v_mat_v_2_d,                             &
                               v_mat_v_2_z
                END INTERFACE  v_mat_v_2
                    INTERFACE v_mat_v_3
              MODULE PROCEDURE v_mat_v_3_d,                             &
                               v_mat_v_3_z
                END INTERFACE  v_mat_v_3
                    INTERFACE v_mat_v_4
              MODULE PROCEDURE v_mat_v_4_d,                             &
                               v_mat_v_4_z
                END INTERFACE  v_mat_v_4
                    INTERFACE v_mat_v_5
              MODULE PROCEDURE v_mat_v_5_d,                             &
                               v_mat_v_5_z
                END INTERFACE  v_mat_v_5
                    INTERFACE v_mat_v_6
              MODULE PROCEDURE v_mat_v_6_d,                             &
                               v_mat_v_6_z
                END INTERFACE  v_mat_v_6
                    INTERFACE v_mat_v_7
              MODULE PROCEDURE v_mat_v_7_d,                             &
                               v_mat_v_7_z
                END INTERFACE  v_mat_v_7
                    INTERFACE v_mat_v_8
              MODULE PROCEDURE v_mat_v_8_d,                             &
                               v_mat_v_8_z
                END INTERFACE  v_mat_v_8
                    INTERFACE v_mat_v_9
              MODULE PROCEDURE v_mat_v_9_d,                             &
                               v_mat_v_9_z
                END INTERFACE  v_mat_v_9
                    INTERFACE v_mat_v_10
              MODULE PROCEDURE v_mat_v_10_d,                            &
                               v_mat_v_10_z
                END INTERFACE  v_mat_v_10
                    INTERFACE v_mat_v_11
              MODULE PROCEDURE v_mat_v_11_d,                            &
                               v_mat_v_11_z
                END INTERFACE  v_mat_v_11
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE v_v_mat_gen
              MODULE PROCEDURE v_v_mat_gen_d,                           &
                               v_v_mat_gen_z
                END INTERFACE  v_v_mat_gen
                    INTERFACE v_v_mat_2
              MODULE PROCEDURE v_v_mat_2_d,                             &
                               v_v_mat_2_z
                END INTERFACE  v_v_mat_2
                    INTERFACE v_v_mat_3
              MODULE PROCEDURE v_v_mat_3_d,                             &
                               v_v_mat_3_z
                END INTERFACE  v_v_mat_3
                    INTERFACE v_v_mat_4
              MODULE PROCEDURE v_v_mat_4_d,                             &
                               v_v_mat_4_z
                END INTERFACE  v_v_mat_4
                    INTERFACE v_v_mat_5
              MODULE PROCEDURE v_v_mat_5_d,                             &
                               v_v_mat_5_z
                END INTERFACE  v_v_mat_5
                    INTERFACE v_v_mat_6
              MODULE PROCEDURE v_v_mat_6_d,                             &
                               v_v_mat_6_z
                END INTERFACE  v_v_mat_6
                    INTERFACE v_v_mat_7
              MODULE PROCEDURE v_v_mat_7_d,                             &
                               v_v_mat_7_z
                END INTERFACE  v_v_mat_7
                    INTERFACE v_v_mat_8
              MODULE PROCEDURE v_v_mat_8_d,                             &
                               v_v_mat_8_z
                END INTERFACE  v_v_mat_8
                    INTERFACE v_v_mat_9
              MODULE PROCEDURE v_v_mat_9_d,                             &
                               v_v_mat_9_z
                END INTERFACE  v_v_mat_9
                    INTERFACE v_v_mat_10
              MODULE PROCEDURE v_v_mat_10_d,                            &
                               v_v_mat_10_z
                END INTERFACE  v_v_mat_10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE v_so_mat_v_gen
              MODULE PROCEDURE v_so_mat_v_gen_d,                        &
                               v_so_mat_v_gen_z
                END INTERFACE  v_so_mat_v_gen
                    INTERFACE v_so_mat_v_2
              MODULE PROCEDURE v_so_mat_v_2_d,                          &
                               v_so_mat_v_2_z
                END INTERFACE  v_so_mat_v_2
                    INTERFACE v_so_mat_v_3
              MODULE PROCEDURE v_so_mat_v_3_d,                          &
                               v_so_mat_v_3_z
                END INTERFACE  v_so_mat_v_3
                    INTERFACE v_so_mat_v_4
              MODULE PROCEDURE v_so_mat_v_4_d,                          &
                               v_so_mat_v_4_z
                END INTERFACE  v_so_mat_v_4
                    INTERFACE v_so_mat_v_5
              MODULE PROCEDURE v_so_mat_v_5_d,                          &
                               v_so_mat_v_5_z
                END INTERFACE  v_so_mat_v_5
                    INTERFACE v_so_mat_v_6
              MODULE PROCEDURE v_so_mat_v_6_d,                          &
                               v_so_mat_v_6_z
                END INTERFACE  v_so_mat_v_6
                    INTERFACE v_so_mat_v_7
              MODULE PROCEDURE v_so_mat_v_7_d,                          &
                               v_so_mat_v_7_z
                END INTERFACE  v_so_mat_v_7
                    INTERFACE v_so_mat_v_8
              MODULE PROCEDURE v_so_mat_v_8_d,                          &
                               v_so_mat_v_8_z
                END INTERFACE  v_so_mat_v_8
                    INTERFACE v_so_mat_v_9
              MODULE PROCEDURE v_so_mat_v_9_d,                          &
                               v_so_mat_v_9_z
                END INTERFACE  v_so_mat_v_9
                    INTERFACE v_so_mat_v_10
              MODULE PROCEDURE v_so_mat_v_10_d,                         &
                               v_so_mat_v_10_z
                END INTERFACE  v_so_mat_v_10
                    INTERFACE v_so_mat_v_11
              MODULE PROCEDURE v_so_mat_v_11_d,                         &
                               v_so_mat_v_11_z
                END INTERFACE  v_so_mat_v_11
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE v_v_so_mat_gen
              MODULE PROCEDURE v_v_so_mat_gen_d,                        &
                               v_v_so_mat_gen_z
                END INTERFACE  v_v_so_mat_gen
                    INTERFACE v_v_so_mat_2
              MODULE PROCEDURE v_v_so_mat_2_d,                          &
                               v_v_so_mat_2_z
                END INTERFACE  v_v_so_mat_2
                    INTERFACE v_v_so_mat_3
              MODULE PROCEDURE v_v_so_mat_3_d,                          &
                               v_v_so_mat_3_z
                END INTERFACE  v_v_so_mat_3
                    INTERFACE v_v_so_mat_4
              MODULE PROCEDURE v_v_so_mat_4_d,                          &
                               v_v_so_mat_4_z
                END INTERFACE  v_v_so_mat_4
                    INTERFACE v_v_so_mat_5
              MODULE PROCEDURE v_v_so_mat_5_d,                          &
                               v_v_so_mat_5_z
                END INTERFACE  v_v_so_mat_5
                    INTERFACE v_v_so_mat_6
              MODULE PROCEDURE v_v_so_mat_6_d,                          &
                               v_v_so_mat_6_z
                END INTERFACE  v_v_so_mat_6
                    INTERFACE v_v_so_mat_7
              MODULE PROCEDURE v_v_so_mat_7_d,                          &
                               v_v_so_mat_7_z
                END INTERFACE  v_v_so_mat_7
                    INTERFACE v_v_so_mat_8
              MODULE PROCEDURE v_v_so_mat_8_d,                          &
                               v_v_so_mat_8_z
                END INTERFACE  v_v_so_mat_8
                    INTERFACE v_v_so_mat_9
              MODULE PROCEDURE v_v_so_mat_9_d,                          &
                               v_v_so_mat_9_z
                END INTERFACE  v_v_so_mat_9
                    INTERFACE v_v_so_mat_10
              MODULE PROCEDURE v_v_so_mat_10_d,                         &
                               v_v_so_mat_10_z
                END INTERFACE  v_v_so_mat_10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   Driving routines for DVR matrix vector multiply
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck finite_element_matrix_multiply
!***begin prologue     finite_element_matrix_multiply     
!***date written       040607   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            matrix multiply routine for DVR matrices.
!
!                     H = H(1,1) + H(2,2) + H(3,3)
!
!     Note that vectors are stored as V(nphy(3),nphy(2),nphy(1)).  Thus matrix
!     vector multiply in two and three dimensions must be done carefully for
!     nphy(2) and nphy(3).  We have( repeated indices summed over), in 3D,
!    
! V(3,2,1) = [ H(1,1) * V(3,2,1) + H(2,2) * V(3,2,1) + H(3,3) * V(3,2,1) ]  
!                                   =
!            [ V(3,2,1) * H(1,1) + V(3,2,1) * H(2,2) + H(3,3) * V(3,2,1) ]
!
! Thus the first mutiply may be done as a matrix V(3*2,1) on H(1,1), 
! the second as V(3,2,1) * H(2,2) with an outer loop over index 1 and 
! the third as a simple matrix vector multiply, H(3,3) * V(3,2*1).
!
!and in 2D,
!
! V(2,1)   = [ H(1,1) * V(2,1) + H(2,2) * V(2,1)  ]  
!                                   =
!            [ V(2,1) * H(1,1) + H(2,2) * V(2,1) ]
!
!
!
!***references
!***routines called    
!***                   
!***                   
!
!***end prologue       finite_element_matrix_multiply
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck finite_element_m_v_d
!***begin prologue     finite_element_m_v_d    
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
!                 V = M * V
!
!------------------------------------------------------------------------------------
!***                   Where M is a DVR matrix.
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
!***end prologue       finite_element_m_v_d                     
!
  SUBROUTINE finite_element_m_v_d(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                 :: nv
  INTEGER                                 :: i
  REAL*8, DIMENSION(n3d,nv)               :: v_in
  REAL*8, DIMENSION(n3d,nv)               :: v_out
!
  v_out=0.d0
!
  IF(spdim == 1) THEN
     CALL dvr_mat_mul(v_in,v_out,nphy(1),nv,1)
     CALL v_diag_mul_1d_d(v_in,v_out,nv)
  ELSE IF(spdim == 2) THEN
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       This code is for the general problem involving nphy(2),
!       that is H(2,2) * V(2,1) 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL dvr_mat_mul(v_in,v_out,nphy(2),nphy(1)*nv,2)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       This code is for the general problem involving nphy(1),
!       V(2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        CALL dvr_mat_mul_2d_d(v_in,v_out,nphy(2),nphy(1),nv,1)
        CALL v_diag_mul_2d_d(v_in,v_out,nv)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ELSE IF(spdim == 3) THEN
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            This code is for the general problem involving nphy(3),
!            H(3,3) * V(3,2,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        CALL dvr_mat_mul(v_in,v_out,nphy(3),nphy(2)*nphy(1)*nv,3)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general case involving the matrix multiply 
!        involving nphy(2), that is, V(3,2,1) * H(2,2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL dvr_mat_mul_2d_d(v_in,v_out,nphy(3),nphy(2),nphy(1)*nv,2)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general problem involving the matrix multiply 
!        for nphy(1), that is, V(3,2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL dvr_mat_mul_2d_d(v_in,v_out,nphy(3)*nphy(2),nphy(1),nv,1)
        CALL v_diag_mul_3d_d(v_in,v_out,nv)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END IF
  END SUBROUTINE finite_element_m_v_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!deck finite_element_m_v_z
!***begin prologue     finite_element_m_v_z    
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
!                 V = M * V
!
!------------------------------------------------------------------------------------
!***                   Where M is a DVR matrix.
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
!***end prologue       finite_element_m_v_z                     
!
  SUBROUTINE finite_element_m_v_z(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                 :: nv
  INTEGER                                 :: i
  COMPLEX*16, DIMENSION(n3d,nv)           :: v_in
  COMPLEX*16, DIMENSION(n3d,nv)           :: v_out
!
  v_out=0.d0
!
  IF(spdim == 1) THEN
     CALL dvr_mat_mul(v_in,v_out,nphy(1),nv,1)
     CALL v_diag_mul_1d_z(v_in,v_out,nv)
  ELSE IF(spdim == 2) THEN
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       This code is for the general problem involving nphy(2),
!       that is H(2,2) * V(2,1) 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL dvr_mat_mul(v_in,v_out,nphy(2),nphy(1)*nv,2)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       This code is for the general problem involving nphy(1),
!       V(2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        CALL dvr_mat_mul_2d_z(v_in,v_out,nphy(2),nphy(1),nv,1)
        CALL v_diag_mul_2d_z(v_in,v_out,nv)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ELSE IF(spdim == 3) THEN
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            This code is for the general problem involving nphy(3),
!            H(3,3) * V(3,2,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        CALL dvr_mat_mul(v_in,v_out,nphy(3),nphy(2)*nphy(1)*nv,3)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general case involving the matrix multiply 
!        involving nphy(2), that is, V(3,2,1) * H(2,2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL dvr_mat_mul_2d_z(v_in,v_out,nphy(3),nphy(2),nphy(1)*nv,2)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general problem involving the matrix multiply 
!        for nphy(1), that is, V(3,2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL dvr_mat_mul_2d_z(v_in,v_out,nphy(3)*nphy(2),nphy(1),nv,1)
        CALL v_diag_mul_3d_z(v_in,v_out,nv)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END IF
  END SUBROUTINE finite_element_m_v_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!deck finite_element_h_v_d
!***begin prologue     finite_element_h_v_d    
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
!                 V = M * V
!
!------------------------------------------------------------------------------------
!***                   Where M is a DVR matrix.
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
!***end prologue       finite_element_h_v_d                     
!
  SUBROUTINE finite_element_h_v_d(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                    :: nv
  INTEGER                                    :: i
  REAL*8, DIMENSION(n3d,nv)                  :: v_in
  REAL*8, DIMENSION(n3d,nv)                  :: v_out
!
  v_out=0.d0
!
  IF(spdim == 1) THEN
     CALL dvr_mat_mul(v_in,v_out,nphy(1),nv,1)
     CALL v_diag_mul_1d_d(v_in,v_out,nv)
  ELSE IF(spdim == 2) THEN
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       This code is for the general problem involving nphy(2),
!       that is H(2,2) * V(2,1) 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL dvr_mat_mul(v_in,v_out,nphy(2),nphy(1)*nv,2)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       This code is for the general problem involving nphy(1),
!       V(2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        CALL dvr_mat_mul_2d_d(v_in,v_out,nphy(2),nphy(1),nv,1)
        CALL v_diag_mul_2d_d(v_in,v_out,nv)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ELSE IF(spdim == 3) THEN
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            This code is for the general problem involving nphy(3),
!            H(3,3) * V(3,2,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        CALL dvr_mat_mul(v_in,v_out,nphy(3),nphy(2)*nphy(1)*nv,3)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general case involving the matrix multiply 
!        involving nphy(2), that is, V(3,2,1) * H(2,2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL dvr_mat_mul_2d_d(v_in,v_out,nphy(3),nphy(2),nphy(1)*nv,2)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general problem involving the matrix multiply 
!        for nphy(1), that is, V(3,2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL dvr_mat_mul_2d_d(v_in,v_out,nphy(3)*nphy(2),nphy(1),nv,1)
        CALL v_diag_mul_3d_d(v_in,v_out,nv)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END IF
  END SUBROUTINE finite_element_h_v_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    Routines to compute the matrix vector
!                    mutiply with DVR Hamiltonian
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!deck dvr_mat_mul_d
!***begin prologue     dvr_mat_mul_d    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine computes the matrix vector multiply 
!***                   needed for a Hamiltonian which is of the FEDVR
!**                    form.  There is no need for a packed matrix here
!**                    as the structure of the H matrix is recognized.
!**                    The double counting of the interface elements must be
!**                    taken care of before calling this routine.  Other
!**                    routines are used to do this by halving the diagonal
!**                    elements to avoid overcounting. 
!------------------------------------------------------------------------------------
!
!                             V = M * V
!
!------------------------------------------------------------------------------------
!***                   Where M is a DVR matrix.
!
!***                   The parameter nj is a dummy and can take on values
!***                   consistent with any dimensional problem.  In a 1D
!***                   case nj=1, in 2D nj=nx and in 3D nj=ny*nx.
!***                   It is assumed that double counting of the diagonals
!***                   has been accounted for by halving the elements at
!***                   the interior, connecting, DVR element points.
!**                    Special routines have been written for matrices
!**                    up to 10*10 for efficiency.
!***references
!***routines called    
!
!***end prologue       dvr_mat_mul_d
!
  SUBROUTINE dvr_mat_mul_d(v,v_scr,ni,nj,index)
  USE dvrprop_global,  ONLY   : mat_reg
  IMPLICIT NONE
  INTEGER                                :: ni, nj, index
  REAL*8, DIMENSION(ni,nj)               :: v
  REAL*8, DIMENSION(ni,nj)               :: v_scr
  INTEGER                                :: i
  INTEGER                                :: locate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  locate = 1
  DO i = 1, num_reg(index)
!
    IF ( nfun_reg(i,index) > 10 ) then
!
!                    General Code
!
         CALL v_mat_v_gen(v(locate:,:),                  &
                          v_scr(locate:,:),              &
                          mat_reg(i,index)%ke_mat_d,     &
                          ni,nj,nfun_reg(i,index))
    ELSE IF (nfun_reg(i,index) == 2) then
! 
!                       Special case for Two by Two
!
        CALL v_mat_v_2(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg(i,index)%ke_mat_d,        &
                       ni,nj)

    ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Special case for Three by Three
!
        CALL v_mat_v_3(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg(i,index)%ke_mat_d,        &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Special case for Four by Four
!
        CALL v_mat_v_4(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg(i,index)%ke_mat_d,        &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Special case for Five by Five
!
        CALL v_mat_v_5(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg(i,index)%ke_mat_d,        &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Special case for Six by Six
!
        CALL v_mat_v_6(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg(i,index)%ke_mat_d,        &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Special case for Seven by Seven
!
        CALL v_mat_v_7(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg(i,index)%ke_mat_d,        &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Special case for Eight by Eight
!
        CALL v_mat_v_8(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg(i,index)%ke_mat_d,        &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Special case for Nine by Nine
!
        CALL v_mat_v_9(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg(i,index)%ke_mat_d,        &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Special case for Ten by Ten
!
        CALL v_mat_v_10(v(locate:,:),                    &
                        v_scr(locate:,:),                &
                        mat_reg(i,index)%ke_mat_d,       &
                        ni,nj)

!    ELSE IF (nfun_reg(i,index) == 11) then
!
!                       Special case for Eleven by Eleven
!
!        CALL v_mat_v_11(v(locate:,:),                    &
!                        v_scr(locate:,:),                &
!                        mat_reg(i,index)%ke_mat_d,       &
!                        ni,nj)
    END IF
    locate = locate + nfun_reg(i,index)  - 1
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE dvr_mat_mul_d
!
!deck dvr_mat_mul_z
!***begin prologue     dvr_mat_mul_z    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine computes the matrix vector multiply 
!***                   needed for a Hamiltonian which is of the FEDVR
!**                    form.  There is no need for a packed matrix here
!**                    as the structure of the M matrix is recognized.
!**                    The double counting of the interface elements must be
!**                    taken care of before calling this routine.  Other
!**                    routines are used to do this by halving the diagonal
!**                    elements to avoid overcounting. 
!------------------------------------------------------------------------------------
!
!                             V = M * V
!
!------------------------------------------------------------------------------------
!***                   Where M is a DVR matrix.
!
!***                   The parameter nj is a dummy and can take on values
!***                   consistent with any dimensional problem.  In a 1D
!***                   case nj=1, in 2D nj=nx and in 3D nj=ny*nx.
!***                   It is assumed that double counting of the diagonals
!***                   has been accounted for by halving the elements at
!***                   the interior, connecting, DVR element points.
!**                    Special routines have been written for matrices
!**                    up to 10*10 for efficiency.
!***references
!***routines called    
!
!***end prologue       dvr_mat_mul                     
!
  SUBROUTINE dvr_mat_mul_z(v,v_scr,ni,nj,index)
  USE dvrprop_global,  ONLY   : mat_reg
  IMPLICIT NONE
  INTEGER                                :: ni, nj, index
  COMPLEX*16, DIMENSION(ni,nj)           :: v
  COMPLEX*16, DIMENSION(ni,nj)           :: v_scr
  INTEGER                                :: i
  INTEGER                                :: locate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  locate = 1
  DO i = 1, num_reg(index)
!
    IF ( nfun_reg(i,index) > 10 ) then
!
!                    General Code
!
         CALL v_mat_v_gen(v(locate:,:),                  &
                          v_scr(locate:,:),              &
                          mat_reg(i,index)%ke_mat_d,     &
                          ni,nj,nfun_reg(i,index))
    ELSE IF (nfun_reg(i,index) == 2) then
! 
!                       Special case for Two by Two
!
        CALL v_mat_v_2(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg(i,index)%ke_mat_d,        &
                       ni,nj)

    ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Special case for Three by Three
!
        CALL v_mat_v_3(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg(i,index)%ke_mat_d,        &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Special case for Four by Four
!
        CALL v_mat_v_4(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg(i,index)%ke_mat_d,        &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Special case for Five by Five
!
        CALL v_mat_v_5(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg(i,index)%ke_mat_d,        &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Special case for Six by Six
!
        CALL v_mat_v_6(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg(i,index)%ke_mat_d,        &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Special case for Seven by Seven
!
        CALL v_mat_v_7(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg(i,index)%ke_mat_d,        &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Special case for Eight by Eight
!
        CALL v_mat_v_8(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg(i,index)%ke_mat_d,        &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Special case for Nine by Nine
!
        CALL v_mat_v_9(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg(i,index)%ke_mat_d,        &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Special case for Ten by Ten
!
        CALL v_mat_v_10(v(locate:,:),                    &
                        v_scr(locate:,:),                &
                        mat_reg(i,index)%ke_mat_d,       &
                        ni,nj)
    END IF
    locate = locate + nfun_reg(i,index)  - 1
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE dvr_mat_mul_z
!
!*deck dvr_mat_mul_2d_d
!***begin prologue     dvr_mat_mul_2d_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine is the driver routine which
!***                   calculates the full two dimensional matrix
!**                    vector multiply using the fact that the matrices
!**                    are separable in each dimension.
!------------------------------------------------------------------------------------
!
!                            V = V * M 
!
!------------------------------------------------------------------------------------
!***                   where M is the regional
!***                   propagator for a FEDVR or FD Hamiltonian. 
!***                   parts.  
!***                   The routine is not needed in 1D problems. 
!***                   In 2D, the ni index runs over the number
!                      of points in the y coordinate while in 3D,
!***                   it runs over the product nz*ny of the number 
!***                   of points in the z and y coordinates.
!
!***references
!***routines called    
!***end prologue       dvr_mat_mul_2d_d
!
  SUBROUTINE dvr_mat_mul_2d_d(v,v_scr,ni,nj,nv,index)
  USE dvrprop_global,  ONLY   : mat_reg
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, nv, index
  REAL*8, DIMENSION(ni,nj,nv)              :: v
  REAL*8, DIMENSION(ni,nj,nv)              :: v_scr
  INTEGER                                  :: i, j
  INTEGER                                  :: locate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  locate = 1
!
  DO i=1, num_reg(index)
     IF ( nfun_reg(i,index) > 10 ) then
!
!                        General Case
!
          DO j=1,nv
             CALL v_v_mat_gen(v(:,locate:,j),                   &
                              v_scr(:,locate:,j),               &
                              mat_reg(i,index)%ke_mat_d,        &
                              ni,nj,nfun_reg(i,index)) 
          END DO
!
     ELSE IF (nfun_reg(i,index) == 2) then
!
!                       Two by Two
! 
          DO j=1,nv
             CALL v_v_mat_2(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg(i,index)%ke_mat_d,          &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Three by Three
!
          DO j=1,nv
             CALL v_v_mat_3(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg(i,index)%ke_mat_d,          &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Four by Four
!
          DO j=1,nv
             CALL v_v_mat_4(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg(i,index)%ke_mat_d,          &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Five by Five
!
          DO j=1,nv
             CALL v_v_mat_5(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg(i,index)%ke_mat_d,          &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Six by Six
!
          DO j=1,nv
             CALL v_v_mat_6(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg(i,index)%ke_mat_d,          &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Seven by Seven
!
          DO j=1,nv
             CALL v_v_mat_7(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg(i,index)%ke_mat_d,          &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Eight by Eight
!
          DO j=1,nv
             CALL v_v_mat_8(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg(i,index)%ke_mat_d,          &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Nine by Nine
!
          DO j=1,nv
             CALL v_v_mat_9(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg(i,index)%ke_mat_d,          &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Ten by Ten
!
          DO j=1,nv
             CALL v_v_mat_10(v(:,locate:,j),                    &
                             v_scr(:,locate:,j),                &
                             mat_reg(i,index)%ke_mat_d,         &
                             ni,nj)
          END DO
     END IF
     locate = locate + nfun_reg(i,index)  - 1
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE dvr_mat_mul_2d_d
!
!*deck dvr_mat_mul_2d_z
!***begin prologue     dvr_mat_mul_2d_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine is the driver routine which
!***                   calculates the full two dimensional matrix
!**                    vector multiply using the fact that the matrices
!**                    are separable in each dimension.
!------------------------------------------------------------------------------------
!
!                            V = V * M 
!
!------------------------------------------------------------------------------------
!***                   where M is the regional
!***                   propagator for a FEDVR or FD Hamiltonian. 
!***                   parts.  
!***                   The routine is not needed in 1D problems. 
!***                   In 2D, the ni index runs over the number
!                      of points in the y coordinate while in 3D,
!***                   it runs over the product nz*ny of the number 
!***                   of points in the z and y coordinates.
!
!***references
!***routines called    
!***end prologue       dvr_mat_mul_2d_z
!
  SUBROUTINE dvr_mat_mul_2d_z(v,v_scr,ni,nj,nv,index)
  USE dvrprop_global,  ONLY   : mat_reg
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, nv, index
  COMPLEX*16, DIMENSION(ni,nj,nv)          :: v
  COMPLEX*16, DIMENSION(ni,nj,nv)          :: v_scr
  INTEGER                                  :: i, j
  INTEGER                                  :: locate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  locate = 1
!
  DO i=1, num_reg(index)
     IF ( nfun_reg(i,index) > 10 ) then
!
!                        General Case
!
          DO j=1,nv
             CALL v_v_mat_gen(v(:,locate:,j),                   &
                              v_scr(:,locate:,j),               &
                              mat_reg(i,index)%ke_mat_d,        &
                              ni,nj,nfun_reg(i,index)) 
          END DO
!
     ELSE IF (nfun_reg(i,index) == 2) then
!
!                       Two by Two
! 
          DO j=1,nv
             CALL v_v_mat_2(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg(i,index)%ke_mat_d,          &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Three by Three
!
          DO j=1,nv
             CALL v_v_mat_3(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg(i,index)%ke_mat_d,          &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Four by Four
!
          DO j=1,nv
             CALL v_v_mat_4(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg(i,index)%ke_mat_d,          &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Five by Five
!
          DO j=1,nv
             CALL v_v_mat_5(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg(i,index)%ke_mat_d,          &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Six by Six
!
          DO j=1,nv
             CALL v_v_mat_6(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg(i,index)%ke_mat_d,          &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Seven by Seven
!
          DO j=1,nv
             CALL v_v_mat_7(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg(i,index)%ke_mat_d,          &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Eight by Eight
!
          DO j=1,nv
             CALL v_v_mat_8(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg(i,index)%ke_mat_d,          &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Nine by Nine
!
          DO j=1,nv
             CALL v_v_mat_9(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg(i,index)%ke_mat_d,          &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Ten by Ten
!
          DO j=1,nv
             CALL v_v_mat_10(v(:,locate:,j),                    &
                             v_scr(:,locate:,j),                &
                             mat_reg(i,index)%ke_mat_d,         &
                             ni,nj)
          END DO
     END IF
     locate = locate + nfun_reg(i,index)  - 1
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE dvr_mat_mul_2d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_diag_mul_1d_d
!***begin prologue     v_diag_mul_1d_d    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description         
!***                   
!***references
!***routines called    
!
!***end prologue       v_diag_mul_1d_d
!
  SUBROUTINE v_diag_mul_1d_d(v_in,v_out,nv)
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: nv
  REAL*8, DIMENSION(nphy(1),nv)          :: v_in
  REAL*8, DIMENSION(nphy(1),nv)          :: v_out
  INTEGER                                :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i = 1, nphy(1)
     v_out(i,:) = v_out(i,:) + grid(1)%v(i) * v_in(i,:)
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END SUBROUTINE v_diag_mul_1d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_diag_mul_2d_d
!***begin prologue     v_diag_mul_2d_d    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description         
!***                   
!***references
!***routines called    
!
!***end prologue       v_diag_mul_2d_d
!
  SUBROUTINE v_diag_mul_2d_d(v_in,v_out,nv)
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: nv
  REAL*8, DIMENSION(nphy(2),nphy(1),nv)  :: v_in
  REAL*8, DIMENSION(nphy(2),nphy(1),nv)  :: v_out
  INTEGER                                :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i = 1, nphy(2)
     v_out(i,:,:) = v_out(i,:,:) + grid(2)%v(i) * v_in(i,:,:)
  END DO
  DO i = 1, nphy(1)
     v_out(:,i,:) = v_out(:,i,:) + grid(1)%v(i) * v_in(:,i,:)
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END SUBROUTINE v_diag_mul_2d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_diag_mul_3d_d
!***begin prologue     v_diag_mul_3d_d    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description         
!***                   
!***references
!***routines called    
!
!***end prologue       v_diag_mul_3d_d
!
  SUBROUTINE v_diag_mul_3d_d(v_in,v_out,nv)
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                        :: nv
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),nv)  :: v_in
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),nv)  :: v_out
  INTEGER                                        :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i = 1, nphy(3)
     v_out(i,:,:,:) = v_out(i,:,:,:) + grid(3)%v(i) * v_in(i,:,:,:)
  END DO
  DO i = 1, nphy(2)
     v_out(:,i,:,:) = v_out(:,i,:,:) + grid(2)%v(i) * v_in(:,i,:,:)
  END DO
  DO i = 1, nphy(1)
     v_out(:,:,i,:) = v_out(:,:,i,:) + grid(1)%v(i) * v_in(:,:,i,:)
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END SUBROUTINE v_diag_mul_3d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_diag_mul_1d_z
!***begin prologue     v_diag_mul_1d_z    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description         
!***                   
!***references
!***routines called    
!
!***end prologue       v_diag_mul_1d_z
!
  SUBROUTINE v_diag_mul_1d_z(v_in,v_out,nv)
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                    :: nv
  COMPLEX*16, DIMENSION(nphy(1),nv)          :: v_in
  COMPLEX*16, DIMENSION(nphy(1),nv)          :: v_out
  INTEGER                                    :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i = 1, nphy(1)
     v_out(i,:) = v_out(i,:) + grid(1)%v(i) * v_in(i,:)
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END SUBROUTINE v_diag_mul_1d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_diag_mul_2d_z
!***begin prologue     v_diag_mul_2d_z    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description         
!***                   
!***references
!***routines called    
!
!***end prologue       v_diag_mul_2d_z
!
  SUBROUTINE v_diag_mul_2d_z(v_in,v_out,nv)
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                    :: nv
  COMPLEX*16, DIMENSION(nphy(2),nphy(1),nv)  :: v_in
  COMPLEX*16, DIMENSION(nphy(2),nphy(1),nv)  :: v_out
  INTEGER                                    :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i = 1, nphy(2)
     v_out(i,:,:) = v_out(i,:,:) + grid(2)%v(i) * v_in(i,:,:)
  END DO
  DO i = 1, nphy(1)
     v_out(:,i,:) = v_out(:,i,:) + grid(1)%v(i) * v_in(:,i,:)
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END SUBROUTINE v_diag_mul_2d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_diag_mul_3d_z
!***begin prologue     v_diag_mul_3d_z    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description         
!***                   
!***references
!***routines called    
!
!***end prologue       v_diag_mul_3d_z
!
  SUBROUTINE v_diag_mul_3d_z(v_in,v_out,nv)
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                            :: nv
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1),nv)  :: v_in
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1),nv)  :: v_out
  INTEGER                                            :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i = 1, nphy(3)
     v_out(i,:,:,:) = v_out(i,:,:,:) + grid(3)%v(i) * v_in(i,:,:,:)
  END DO
  DO i = 1, nphy(2)
     v_out(:,i,:,:) = v_out(:,i,:,:) + grid(2)%v(i) * v_in(:,i,:,:)
  END DO
  DO i = 1, nphy(1)
     v_out(:,:,i,:) = v_out(:,:,i,:) + grid(1)%v(i) * v_in(:,:,i,:)
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END SUBROUTINE v_diag_mul_3d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  The routines that do the real work
!                  of the matrix vector mutiply 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**begin prologue     v_mat_v
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       v_mat_v
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_gen_d
!***begin prologue     v_mat_v_gen_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propator matrix vector multiplies
!***                   for a general FEDVR Hamiltonian.
!***                   The parameter nj is a dummy and can take on values
!***                   consistent with any dimensional problem.  In a 1D
!***                   case nj=1, in 2D nj=nx and in 3D nj=ny*nx.
!
!***references
!***routines called    ebcxx, ambcxx, apbcxx
!***end prologue       
!
  SUBROUTINE v_mat_v_gen_d(v,             &
                           v_scr,         &
                           dvr_mat,       &
                           ni,nj,nk)
  IMPLICIT NONE
  INTEGER                              :: ni, nj, nk
  INTEGER                              :: i, k
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(nk,nk)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   DO i=1,nk
      DO k=1,nk
         v_scr(k,:) = v_scr(k,:) + dvr_mat(k,i) * v(i,:) 
      END DO
   END DO
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE v_mat_v_gen_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_2_d
!***begin prologue     v_mat_v_2_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies for
!***                   the special case of 2*2 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_2_d
!
  SUBROUTINE v_mat_v_2_d(v,         &
                         v_scr,     &
                         dvr_mat,   &
                         ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(2,2)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = v_scr(1,:)              &
                            +          &
               dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   

  v_scr(2,:) = v_scr(2,:)              &
                            +          &
               dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_2_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_3_d
!***begin prologue     v_mat_v_3_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_3_d
!
  SUBROUTINE v_mat_v_3_d(v,            &
                         v_scr,        &
                         dvr_mat,      &
                         ni,nj)
  IMPLICIT NONE
  INTEGER                                :: ni, nj
  REAL*8, DIMENSION(ni,nj)               :: v
  REAL*8, DIMENSION(ni,nj)               :: v_scr
  REAL*8, DIMENSION(3,3)                 :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = v_scr(1,:)              &
                            +          &
               dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   

  v_scr(2,:) = v_scr(2,:)              &
                            +          &
               dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   

  v_scr(3,:) = v_scr(3,:)              &
                            +          &
               dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_3_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_4_d
!***begin prologue     v_mat_v_4_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 4*4 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_4_d
!
  SUBROUTINE v_mat_v_4_d(v,             &
                         v_scr,         &
                         dvr_mat,       &
                         ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(4,4)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = v_scr(1,:)              &
                            +          &
               dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   &
                            +          &   
               dvr_mat(1,4) * v(4,:)   

  v_scr(2,:) = v_scr(2,:)              &
                            +          &
               dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   &
                            +          &
               dvr_mat(2,4) * v(4,:)   

  v_scr(3,:) = v_scr(3,:)              &
                            +          &
               dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   &
                            +          &
               dvr_mat(3,4) * v(4,:)   

  v_scr(4,:) = v_scr(4,:)              &
                            +          &
               dvr_mat(4,1) * v(1,:)   &
                            +          &
               dvr_mat(4,2) * v(2,:)   &
                            +          &
               dvr_mat(4,3) * v(3,:)   &
                            +          &
               dvr_mat(4,4) * v(4,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_4_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_5_d
!***begin prologue     v_mat_v_5_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_5_d
!
  SUBROUTINE v_mat_v_5_d(v,          &
                         v_scr,      &
                         dvr_mat,    &
                         ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(5,5)               :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = v_scr(1,:)              &
                            +          &
               dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   &
                            +          &   
               dvr_mat(1,4) * v(4,:)   &
                            +          &
               dvr_mat(1,5) * v(5,:)    

  v_scr(2,:) = v_scr(2,:)              &
                            +          &
               dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   &
                            +          &
               dvr_mat(2,4) * v(4,:)   &
                            +          &
               dvr_mat(2,5) * v(5,:)   

  v_scr(3,:) = v_scr(3,:)              &
                            +          &
               dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   &
                            +          &
               dvr_mat(3,4) * v(4,:)   &
                            +          &
               dvr_mat(3,5) * v(5,:)   

  v_scr(4,:) = v_scr(4,:)              &
                            +          &
               dvr_mat(4,1) * v(1,:)   &
                            +          &
               dvr_mat(4,2) * v(2,:)   &
                            +          &
               dvr_mat(4,3) * v(3,:)   &
                            +          &
               dvr_mat(4,4) * v(4,:)   &
                            +          &
               dvr_mat(4,5) * v(5,:)    
  v_scr(5,:) = v_scr(5,:)              &
                            +          &
               dvr_mat(5,1) * v(1,:)   &
                            +          &
               dvr_mat(5,2) * v(2,:)   &
                            +          &
               dvr_mat(5,3) * v(3,:)   &
                            +          &
               dvr_mat(5,4) * v(4,:)   &
                            +          &
               dvr_mat(5,5) * v(5,:)    
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_5_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_6_d
!***begin prologue     v_mat_v_6_d
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_6_d
!
  SUBROUTINE v_mat_v_6_d(v,           &    
                         v_scr,       &
                         dvr_mat,     &
                         ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(6,6)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = v_scr(1,:)              &
                            +          &
               dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   &
                            +          &   
               dvr_mat(1,4) * v(4,:)   &
                            +          &
               dvr_mat(1,5) * v(5,:)   &   
                            +          &
               dvr_mat(1,6) * v(6,:)    

  v_scr(2,:) = v_scr(2,:)              &
                            +          &   
               dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   &
                            +          &   
               dvr_mat(2,4) * v(4,:)   &
                            +          &
               dvr_mat(2,5) * v(5,:)   &   
                            +          &
               dvr_mat(2,6) * v(6,:)      

  v_scr(3,:) = v_scr(3,:)              &
                            +          & 
               dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   &
                            +          &   
               dvr_mat(3,4) * v(4,:)   &
                            +          &
               dvr_mat(3,5) * v(5,:)   &   
                            +          &
               dvr_mat(3,6) * v(6,:)      

  v_scr(4,:) = v_scr(4,:)              &
                            +          &  
               dvr_mat(4,1) * v(1,:)   &
                            +          &
               dvr_mat(4,2) * v(2,:)   &
                            +          &
               dvr_mat(4,3) * v(3,:)   &
                            +          &   
               dvr_mat(4,4) * v(4,:)   &
                            +          &
               dvr_mat(4,5) * v(5,:)   &   
                            +          &
               dvr_mat(4,6) * v(6,:)      

  v_scr(5,:) = v_scr(5,:)              &
                            +          &  
               dvr_mat(5,1) * v(1,:)   &
                            +          &
               dvr_mat(5,2) * v(2,:)   &
                            +          &
               dvr_mat(5,3) * v(3,:)   &
                            +          &   
               dvr_mat(5,4) * v(4,:)   &
                            +          &
               dvr_mat(5,5) * v(5,:)   &   
                            +          &
               dvr_mat(5,6) * v(6,:)      

  v_scr(6,:) = v_scr(6,:)              &
                            +          &  
               dvr_mat(6,1) * v(1,:)   &
                            +          &
               dvr_mat(6,2) * v(2,:)   &
                            +          &
               dvr_mat(6,3) * v(3,:)   &
                            +          &   
               dvr_mat(6,4) * v(4,:)   &
                            +          &
               dvr_mat(6,5) * v(5,:)   &   
                            +          &
               dvr_mat(6,6) * v(6,:)      
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_6_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_7_d
!***begin prologue     v_mat_v_7_d  
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_7_d
!
  SUBROUTINE v_mat_v_7_d(v,           &    
                         v_scr,       &
                         dvr_mat,     &
                         ni,nj)
                     
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(7,7)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = v_scr(1,:)              &
                            +          &
               dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   &
                            +          &   
               dvr_mat(1,4) * v(4,:)   &
                            +          &
               dvr_mat(1,5) * v(5,:)   &   
                            +          &
               dvr_mat(1,6) * v(6,:)   &   
                            +          &
               dvr_mat(1,7) * v(7,:)      

  v_scr(2,:) = v_scr(2,:)              &
                            +          &
               dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   &
                            +          &   
               dvr_mat(2,4) * v(4,:)   &
                            +          &
               dvr_mat(2,5) * v(5,:)   &   
                            +          &
               dvr_mat(2,6) * v(6,:)   &   
                            +          &
               dvr_mat(2,7) * v(7,:)    

  v_scr(3,:) = v_scr(3,:)              &
                            +          &  
               dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   &
                            +          &   
               dvr_mat(3,4) * v(4,:)   &
                            +          &
               dvr_mat(3,5) * v(5,:)   &   
                            +          &
               dvr_mat(3,6) * v(6,:)   &   
                            +          &
               dvr_mat(3,7) * v(7,:)    

  v_scr(4,:) = v_scr(4,:)              &
                            +          &  
               dvr_mat(4,1) * v(1,:)   &
                            +          &
               dvr_mat(4,2) * v(2,:)   &
                            +          &
               dvr_mat(4,3) * v(3,:)   &
                            +          &   
               dvr_mat(4,4) * v(4,:)   &
                            +          &
               dvr_mat(4,5) * v(5,:)   &   
                            +          &
               dvr_mat(4,6) * v(6,:)   &   
                            +          &
                dvr_mat(4,7) * v(7,:)   

  v_scr(5,:) = v_scr(5,:)              &
                            +          &   
               dvr_mat(5,1) * v(1,:)   &
                            +          &
               dvr_mat(5,2) * v(2,:)   &
                            +          &
               dvr_mat(5,3) * v(3,:)   &
                            +          &   
               dvr_mat(5,4) * v(4,:)   &
                            +          &
               dvr_mat(5,5) * v(5,:)   &   
                            +          &
               dvr_mat(5,6) * v(6,:)   &   
                            +          &
               dvr_mat(5,7) * v(7,:)    

  v_scr(6,:) = v_scr(6,:)              &
                            +          &  
               dvr_mat(6,1) * v(1,:)   &
                            +          &
               dvr_mat(6,2) * v(2,:)   &
                            +          &
               dvr_mat(6,3) * v(3,:)   &
                            +          &   
               dvr_mat(6,4) * v(4,:)   &
                            +          &
               dvr_mat(6,5) * v(5,:)   &   
                            +          &
               dvr_mat(6,6) * v(6,:)   &   
                            +          &
               dvr_mat(6,7) * v(7,:)    

  v_scr(7,:) = v_scr(7,:)              &
                            +          &  
               dvr_mat(7,1) * v(1,:)   &
                            +          &
               dvr_mat(7,2) * v(2,:)   &
                            +          &
               dvr_mat(7,3) * v(3,:)   &
                            +          &   
               dvr_mat(7,4) * v(4,:)   &
                            +          &
               dvr_mat(7,5) * v(5,:)   &   
                            +          &
               dvr_mat(7,6) * v(6,:)   &   
                            +          &
               dvr_mat(7,7) * v(7,:)      
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_7_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_8_d
!***begin prologue     v_mat_v_8_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_8_d
!
  SUBROUTINE v_mat_v_8_d(v,             &
                         v_scr,         &
                         dvr_mat,       &
                         ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(8,8)               :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  v_scr(1,:) = v_scr(1,:)                 &
                            +             &
               dvr_mat(1,1) * v(1,:)      &           
                            +             &
               dvr_mat(1,2) * v(2,:)      &
                            +             &
               dvr_mat(1,3) * v(3,:)      &
                            +             &   
               dvr_mat(1,4) * v(4,:)      &
                            +             &
               dvr_mat(1,5) * v(5,:)      &   
                            +             &
               dvr_mat(1,6) * v(6,:)      &   
                            +             &
               dvr_mat(1,7) * v(7,:)      &   
                            +             &
               dvr_mat(1,8) * v(8,:)      

  v_scr(2,:) = v_scr(2,:)                 &
                               +          &
               dvr_mat(2,1) * v(1,:)      &
                            +             &
               dvr_mat(2,2) * v(2,:)      &
                            +             &
               dvr_mat(2,3) * v(3,:)      &
                            +             &   
               dvr_mat(2,4) * v(4,:)      &
                            +             &
               dvr_mat(2,5) * v(5,:)      &   
                            +             &
               dvr_mat(2,6) * v(6,:)      &   
                            +             &
               dvr_mat(2,7) * v(7,:)      &   
                            +             &
               dvr_mat(2,8) * v(8,:) 

  v_scr(3,:) = v_scr(3,:)                 &
                            +             &     
               dvr_mat(3,1) * v(1,:)      &
                            +             &
               dvr_mat(3,2) * v(2,:)      &
                            +             &
               dvr_mat(3,3) * v(3,:)      &
                            +             &   
               dvr_mat(3,4) * v(4,:)      &
                            +             &
               dvr_mat(3,5) * v(5,:)      &   
                            +             &
               dvr_mat(3,6) * v(6,:)      &   
                            +             &
               dvr_mat(3,7) * v(7,:)      &   
                            +             &
               dvr_mat(3,8) * v(8,:)      

  v_scr(4,:) = v_scr(4,:)                 &
                            +             &
               dvr_mat(4,1) * v(1,:)      &
                            +             &
               dvr_mat(4,2) * v(2,:)      &
                            +             &
               dvr_mat(4,3) * v(3,:)      &
                            +             &   
               dvr_mat(4,4) * v(4,:)      &
                            +             &
               dvr_mat(4,5) * v(5,:)      &   
                            +             &
               dvr_mat(4,6) * v(6,:)      &   
                            +             &
               dvr_mat(4,7) * v(7,:)      &   
                            +             &
               dvr_mat(4,8) * v(8,:)      

  v_scr(5,:) = v_scr(5,:)                &
                            +             &
               dvr_mat(5,1) * v(1,:)      &
                            +             &
               dvr_mat(5,2) * v(2,:)      &
                            +             &
               dvr_mat(5,3) * v(3,:)      &
                            +             &   
               dvr_mat(5,4) * v(4,:)      &
                            +             &
               dvr_mat(5,5) * v(5,:)      &   
                            +             &
               dvr_mat(5,6) * v(6,:)      &   
                            +             &
               dvr_mat(5,7) * v(7,:)      &   
                            +             &
               dvr_mat(5,8) * v(8,:)     
   
  v_scr(6,:) = v_scr(6,:)                 &
                            +             &
               dvr_mat(6,1) * v(1,:)      &
                            +             &
               dvr_mat(6,2) * v(2,:)      &
                            +             &
               dvr_mat(6,3) * v(3,:)      &
                            +             &   
               dvr_mat(6,4) * v(4,:)      &
                            +             &
               dvr_mat(6,5) * v(5,:)      &   
                            +             &
               dvr_mat(6,6) * v(6,:)      &   
                            +             &
               dvr_mat(6,7) * v(7,:)      &   
                            +             &
               dvr_mat(6,8) * v(8,:)        

  v_scr(7,:) = v_scr(7,:)                 &
                            +             &
               dvr_mat(7,1) * v(1,:)      &
                            +             &
               dvr_mat(7,2) * v(2,:)      &
                            +             &
               dvr_mat(7,3) * v(3,:)      &
                            +             &   
               dvr_mat(7,4) * v(4,:)      &
                            +             &
               dvr_mat(7,5) * v(5,:)      &   
                            +             &
               dvr_mat(7,6) * v(6,:)      &   
                            +             &
               dvr_mat(7,7) * v(7,:)      &   
                            +             &
               dvr_mat(7,8) * v(8,:)        

  v_scr(8,:) = v_scr(8,:)                 &
                            +             &
               dvr_mat(8,1) * v(1,:)      &
                            +             &
               dvr_mat(8,2) * v(2,:)      &
                            +             &
               dvr_mat(8,3) * v(3,:)      &
                            +             &   
               dvr_mat(8,4) * v(4,:)      &
                            +             &
               dvr_mat(8,5) * v(5,:)      &  
                            +             &
               dvr_mat(8,6) * v(6,:)      &   
                            +             &
               dvr_mat(8,7) * v(7,:)      &   
                            +             &
               dvr_mat(8,8) * v(8,:)        
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_8_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_9_d
!***begin prologue     v_mat_v_9_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_9_d
!
  SUBROUTINE v_mat_v_9_d(v,           &
                         v_scr,       &
                         dvr_mat,     &
                         ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(9,9)               :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = v_scr(1,:)                 &
                            +             &
               dvr_mat(1,1) * v(1,:)      &
                            +             &
               dvr_mat(1,2) * v(2,:)      &
                            +             &
               dvr_mat(1,3) * v(3,:)      &
                            +             &   
               dvr_mat(1,4) * v(4,:)      &
                            +             &
               dvr_mat(1,5) * v(5,:)      &   
                            +             &
               dvr_mat(1,6) * v(6,:)      &   
                            +             &
               dvr_mat(1,7) * v(7,:)      &   
                            +             &
               dvr_mat(1,8) * v(8,:)      &   
                            +             &
               dvr_mat(1,9) * v(9,:)    

  v_scr(2,:) = v_scr(2,:)                 &
                            +             &  
               dvr_mat(2,1) * v(1,:)      &
                            +             &
               dvr_mat(2,2) * v(2,:)      &
                            +             &
               dvr_mat(2,3) * v(3,:)      &
                            +             &   
               dvr_mat(2,4) * v(4,:)      &
                            +             &
               dvr_mat(2,5) * v(5,:)      &
                            +             &
               dvr_mat(2,6) * v(6,:)      &   
                            +             &
               dvr_mat(2,7) * v(7,:)      &   
                            +             &
               dvr_mat(2,8) * v(8,:)      &   
                            +             &
               dvr_mat(2,9) * v(9,:)    

  v_scr(3,:) = v_scr(3,:)                 &
                            +             &  
               dvr_mat(3,1) * v(1,:)      &
                            +             &
               dvr_mat(3,2) * v(2,:)      &
                            +             &
               dvr_mat(3,3) * v(3,:)      &
                            +             &   
               dvr_mat(3,4) * v(4,:)      &
                            +             &
               dvr_mat(3,5) * v(5,:)      &   
                            +             &
               dvr_mat(3,6) * v(6,:)      &   
                            +             &
               dvr_mat(3,7) * v(7,:)      &   
                            +             &
               dvr_mat(3,8) * v(8,:)      &   
                            +             &
               dvr_mat(3,9) * v(9,:)      

  v_scr(4,:) = v_scr(4,:)                 &
                            +             &
               dvr_mat(4,1) * v(1,:)      &
                            +             &
               dvr_mat(4,2) * v(2,:)      &
                            +             &
               dvr_mat(4,3) * v(3,:)      &
                            +             &   
               dvr_mat(4,4) * v(4,:)      &
                            +             &
               dvr_mat(4,5) * v(5,:)      &   
                            +             &
               dvr_mat(4,6) * v(6,:)      &   
                            +             &
               dvr_mat(4,7) * v(7,:)      &   
                            +             &
               dvr_mat(4,8) * v(8,:)      &   
                            +             &
               dvr_mat(4,9) * v(9,:)      

  v_scr(5,:) = v_scr(5,:)                 &
                            +             &
               dvr_mat(5,1) * v(1,:)      &
                            +             &
               dvr_mat(5,2) * v(2,:)      &
                            +             &
               dvr_mat(5,3) * v(3,:)      &
                            +             &   
               dvr_mat(5,4) * v(4,:)      &
                            +             &
               dvr_mat(5,5) * v(5,:)      &   
                            +             &
               dvr_mat(5,6) * v(6,:)      &   
                            +             &
               dvr_mat(5,7) * v(7,:)      &   
                            +             &
               dvr_mat(5,8) * v(8,:)      &   
                            +             &
               dvr_mat(5,9) * v(9,:)      

  v_scr(6,:) = v_scr(6,:)                 &
                            +             &
               dvr_mat(6,1) * v(1,:)      &
                            +             &
               dvr_mat(6,2) * v(2,:)      &
                            +             &
               dvr_mat(6,3) * v(3,:)      &
                            +             &   
               dvr_mat(6,4) * v(4,:)      &
                            +             &
               dvr_mat(6,5) * v(5,:)      &   
                            +             &
               dvr_mat(6,6) * v(6,:)      &   
                            +             &
               dvr_mat(6,7) * v(7,:)      &   
                            +             &
               dvr_mat(6,8) * v(8,:)      &   
                            +             &
               dvr_mat(6,9) * v(9,:)      

  v_scr(7,:) = v_scr(7,:)                 &
                            +             &
               dvr_mat(7,1) * v(1,:)      &
                            +             &
               dvr_mat(7,2) * v(2,:)      &
                            +             &
               dvr_mat(7,3) * v(3,:)      &
                            +             &   
               dvr_mat(7,4) * v(4,:)      &
                            +             &
               dvr_mat(7,5) * v(5,:)      &   
                            +             &
               dvr_mat(7,6) * v(6,:)      &  
                            +             &
               dvr_mat(7,7) * v(7,:)      &   
                            +             &
               dvr_mat(7,8) * v(8,:)      &   
                            +             &
               dvr_mat(7,9) * v(9,:)      

  v_scr(8,:) = v_scr(8,:)                 &
                            +             &
               dvr_mat(8,1) * v(1,:)      &
                            +             &
               dvr_mat(8,2) * v(2,:)      &
                            +             &
               dvr_mat(8,3) * v(3,:)      &
                            +             &   
               dvr_mat(8,4) * v(4,:)      &
                            +             &
               dvr_mat(8,5) * v(5,:)      &   
                            +             &
               dvr_mat(8,6) * v(6,:)      &   
                            +             &
               dvr_mat(8,7) * v(7,:)      &   
                            +             &
               dvr_mat(8,8) * v(8,:)      &   
                            +             &
               dvr_mat(8,9) * v(9,:)      

  v_scr(9,:) = v_scr(9,:)                 &
                            +             &
               dvr_mat(9,1) * v(1,:)      &
                            +             &
               dvr_mat(9,2) * v(2,:)      &
                            +             &
               dvr_mat(9,3) * v(3,:)      &
                            +             &   
               dvr_mat(9,4) * v(4,:)      &
                            +             &
               dvr_mat(9,5) * v(5,:)      &   
                            +             &
               dvr_mat(9,6) * v(6,:)      &   
                            +             &
               dvr_mat(9,7) * v(7,:)      &   
                            +             &
               dvr_mat(9,8) * v(8,:)      &   
                            +             &
               dvr_mat(9,9) * v(9,:)      
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_9_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_10_d
!***begin prologue     v_mat_v_10_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_10_d
!
  SUBROUTINE v_mat_v_10_d(v,          &
                          v_scr,      &
                          dvr_mat,    &
                          ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(10,10)             :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = v_scr(1,:)                 &
                            +             &
               dvr_mat(1,1) * v(1,:)      &
                            +             &
               dvr_mat(1,2) * v(2,:)      &
                            +             &
               dvr_mat(1,3) * v(3,:)      &
                            +             &   
               dvr_mat(1,4) * v(4,:)      &
                            +             &
               dvr_mat(1,5) * v(5,:)      &   
                            +             &
               dvr_mat(1,6) * v(6,:)      &   
                            +             &
               dvr_mat(1,7) * v(7,:)      &   
                            +             &
               dvr_mat(1,8) * v(8,:)      &   
                            +             &
               dvr_mat(1,9) * v(9,:)      &   
                            +             &
               dvr_mat(1,10) * v(10,:)    

  v_scr(2,:) = v_scr(2,:)                 &
                            +             &
               dvr_mat(2,1) * v(1,:)      &
                            +             &
               dvr_mat(2,2) * v(2,:)      &
                            +             &
               dvr_mat(2,3) * v(3,:)      &
                            +             &   
               dvr_mat(2,4) * v(4,:)      &
                            +             &
               dvr_mat(2,5) * v(5,:)      &   
                            +             &
               dvr_mat(2,6) * v(6,:)      &   
                            +             &
               dvr_mat(2,7) * v(7,:)      &   
                            +             &
               dvr_mat(2,8) * v(8,:)      &   
                            +             &
               dvr_mat(2,9) * v(9,:)      &   
                            +             &
               dvr_mat(2,10) * v(10,:)    

  v_scr(3,:) = v_scr(3,:)                 &
                            +             &
               dvr_mat(3,1) * v(1,:)      &
                            +             &
               dvr_mat(3,2) * v(2,:)      &
                            +             &
               dvr_mat(3,3) * v(3,:)      &
                            +             &   
               dvr_mat(3,4) * v(4,:)      &
                            +             &
               dvr_mat(3,5) * v(5,:)      &   
                            +             &
               dvr_mat(3,6) * v(6,:)      &   
                            +             &
               dvr_mat(3,7) * v(7,:)      &   
                            +             &
               dvr_mat(3,8) * v(8,:)      &   
                            +             &
               dvr_mat(3,9) * v(9,:)      &   
                            +             &
               dvr_mat(3,10) * v(10,:)    

  v_scr(4,:) = v_scr(4,:)                 &
                            +             &
               dvr_mat(4,1) * v(1,:)      &
                            +             &
               dvr_mat(4,2) * v(2,:)      &
                            +             &
               dvr_mat(4,3) * v(3,:)      &
                            +             &   
               dvr_mat(4,4) * v(4,:)      &
                            +             &
               dvr_mat(4,5) * v(5,:)      &   
                            +             &
               dvr_mat(4,6) * v(6,:)      &   
                            +             &
               dvr_mat(4,7) * v(7,:)      &   
                            +             &
               dvr_mat(4,8) * v(8,:)      &   
                            +             &
               dvr_mat(4,9) * v(9,:)      &   
                            +             &
               dvr_mat(4,10) * v(10,:)    

  v_scr(5,:) = v_scr(5,:)                 &
                            +             &
               dvr_mat(5,1) * v(1,:)      &
                            +             &
               dvr_mat(5,2) * v(2,:)      &
                            +             &
               dvr_mat(5,3) * v(3,:)      &
                            +             &   
               dvr_mat(5,4) * v(4,:)      &
                            +             &
               dvr_mat(5,5) * v(5,:)      &   
                            +             &
               dvr_mat(5,6) * v(6,:)      &   
                            +             &
               dvr_mat(5,7) * v(7,:)      &   
                            +             &
               dvr_mat(5,8) * v(8,:)      &   
                            +             &
               dvr_mat(5,9) * v(9,:)      &   
                            +             &
               dvr_mat(5,10) * v(10,:)  

  v_scr(6,:) = v_scr(6,:)                 &
                            +             &  
               dvr_mat(6,1) * v(1,:)      &
                            +             &
               dvr_mat(6,2) * v(2,:)      &
                            +             &
               dvr_mat(6,3) * v(3,:)      &
                            +             &   
               dvr_mat(6,4) * v(4,:)      &
                            +             &
               dvr_mat(6,5) * v(5,:)      &   
                            +             &
               dvr_mat(6,6) * v(6,:)      &   
                            +             &
               dvr_mat(6,7) * v(7,:)      &   
                            +             &
               dvr_mat(6,8) * v(8,:)      &   
                            +             &
               dvr_mat(6,9) * v(9,:)      &  
                            +             &
               dvr_mat(6,10) * v(10,:)     

  v_scr(7,:) = v_scr(7,:)                 &
                            +             &
               dvr_mat(7,1) * v(1,:)      &
                            +             &
               dvr_mat(7,2) * v(2,:)      &
                            +             &
               dvr_mat(7,3) * v(3,:)      &
                            +             &   
               dvr_mat(7,4) * v(4,:)      &
                            +             &
               dvr_mat(7,5) * v(5,:)      &   
                            +             &
               dvr_mat(7,6) * v(6,:)      &   
                            +             &
               dvr_mat(7,7) * v(7,:)      &   
                            +             &
               dvr_mat(7,8) * v(8,:)      &   
                            +             &
               dvr_mat(7,9) * v(9,:)      &   
                            +             &
               dvr_mat(7,10) * v(10,:)   

  v_scr(8,:) = v_scr(8,:)                 &
                            +             &
               dvr_mat(8,1) * v(1,:)      &
                            +             &
               dvr_mat(8,2) * v(2,:)      &
                            +             &
               dvr_mat(8,3) * v(3,:)      &
                            +             &   
               dvr_mat(8,4) * v(4,:)      &
                            +             &
               dvr_mat(8,5) * v(5,:)      &  
                            +             &
               dvr_mat(8,6) * v(6,:)      &   
                            +             &
               dvr_mat(8,7) * v(7,:)      &   
                            +             &
               dvr_mat(8,8) * v(8,:)      &   
                            +             &
               dvr_mat(8,9) * v(9,:)      &   
                            +             &
               dvr_mat(8,10) * v(10,:)  

  v_scr(9,:) = v_scr(9,:)                 &
                            +             &    
               dvr_mat(9,1) * v(1,:)      &
                            +             &
               dvr_mat(9,2) * v(2,:)      &
                            +             &
               dvr_mat(9,3) * v(3,:)      &
                            +             &   
               dvr_mat(9,4) * v(4,:)      &
                            +             &
               dvr_mat(9,5) * v(5,:)      &   
                            +             &
               dvr_mat(9,6) * v(6,:)      &   
                            +             &
               dvr_mat(9,7) * v(7,:)      &   
                            +             &
               dvr_mat(9,8) * v(8,:)      &   
                            +             &
               dvr_mat(9,9) * v(9,:)      &   
                            +             &
               dvr_mat(9,10) * v(10,:)    

  v_scr(10,:) = v_scr(10,:)              &
                            +            &
               dvr_mat(10,1) * v(1,:)    &
                              +          &
               dvr_mat(10,2) * v(2,:)    &
                              +          &
               dvr_mat(10,3) * v(3,:)    &
                              +          &   
               dvr_mat(10,4) * v(4,:)    &
                              +          &
               dvr_mat(10,5) * v(5,:)    &   
                              +          &
               dvr_mat(10,6) * v(6,:)    &   
                              +          &
               dvr_mat(10,7) * v(7,:)    &   
                              +          &
               dvr_mat(10,8) * v(8,:)    &   
                              +          &
               dvr_mat(10,9) * v(9,:)    &   
                              +          &
               dvr_mat(10,10) * v(10,:)    
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_10_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_11_d
!***begin prologue     v_mat_v_11_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_11_d
!
  SUBROUTINE v_mat_v_11_d(v,          &
                          v_scr,      &
                          dvr_mat,    &
                          ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(11,11)             :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = v_scr(1,:)                 &
                            +             &
               dvr_mat(1,1) * v(1,:)      &
                            +             &
               dvr_mat(1,2) * v(2,:)      &
                            +             &
               dvr_mat(1,3) * v(3,:)      &
                            +             &   
               dvr_mat(1,4) * v(4,:)      &
                            +             &
               dvr_mat(1,5) * v(5,:)      &   
                            +             &
               dvr_mat(1,6) * v(6,:)      &   
                            +             &
               dvr_mat(1,7) * v(7,:)      &   
                            +             &
               dvr_mat(1,8) * v(8,:)      &   
                            +             &
               dvr_mat(1,9) * v(9,:)      &   
                            +             &
               dvr_mat(1,10) * v(10,:)    &
                            +             &
               dvr_mat(1,11) * v(11,:)    

  v_scr(2,:) = v_scr(2,:)                 &
                            +             &
               dvr_mat(2,1) * v(1,:)      &
                            +             &
               dvr_mat(2,2) * v(2,:)      &
                            +             &
               dvr_mat(2,3) * v(3,:)      &
                            +             &   
               dvr_mat(2,4) * v(4,:)      &
                            +             &
               dvr_mat(2,5) * v(5,:)      &   
                            +             &
               dvr_mat(2,6) * v(6,:)      &   
                            +             &
               dvr_mat(2,7) * v(7,:)      &   
                            +             &
               dvr_mat(2,8) * v(8,:)      &   
                            +             &
               dvr_mat(2,9) * v(9,:)      &   
                            +             &
               dvr_mat(2,10) * v(10,:)    &
                            +             &
               dvr_mat(2,11) * v(11,:)    

  v_scr(3,:) = v_scr(3,:)                 &
                            +             &
               dvr_mat(3,1) * v(1,:)      &
                            +             &
               dvr_mat(3,2) * v(2,:)      &
                            +             &
               dvr_mat(3,3) * v(3,:)      &
                            +             &   
               dvr_mat(3,4) * v(4,:)      &
                            +             &
               dvr_mat(3,5) * v(5,:)      &   
                            +             &
               dvr_mat(3,6) * v(6,:)      &   
                            +             &
               dvr_mat(3,7) * v(7,:)      &   
                            +             &
               dvr_mat(3,8) * v(8,:)      &   
                            +             &
               dvr_mat(3,9) * v(9,:)      &   
                            +             &
               dvr_mat(3,10) * v(10,:)    &
                            +             &
               dvr_mat(3,11) * v(11,:)    

  v_scr(4,:) = v_scr(4,:)                 &
                            +             &
               dvr_mat(4,1) * v(1,:)      &
                            +             &
               dvr_mat(4,2) * v(2,:)      &
                            +             &
               dvr_mat(4,3) * v(3,:)      &
                            +             &   
               dvr_mat(4,4) * v(4,:)      &
                            +             &
               dvr_mat(4,5) * v(5,:)      &   
                            +             &
               dvr_mat(4,6) * v(6,:)      &   
                            +             &
               dvr_mat(4,7) * v(7,:)      &   
                            +             &
               dvr_mat(4,8) * v(8,:)      &   
                            +             &
               dvr_mat(4,9) * v(9,:)      &   
                            +             &
               dvr_mat(4,10) * v(10,:)    &
                            +             &
               dvr_mat(4,11) * v(11,:)    

  v_scr(5,:) = v_scr(5,:)                 &
                            +             &
               dvr_mat(5,1) * v(1,:)      &
                            +             &
               dvr_mat(5,2) * v(2,:)      &
                            +             &
               dvr_mat(5,3) * v(3,:)      &
                            +             &   
               dvr_mat(5,4) * v(4,:)      &
                            +             &
               dvr_mat(5,5) * v(5,:)      &   
                            +             &
               dvr_mat(5,6) * v(6,:)      &   
                            +             &
               dvr_mat(5,7) * v(7,:)      &   
                            +             &
               dvr_mat(5,8) * v(8,:)      &   
                            +             &
               dvr_mat(5,9) * v(9,:)      &   
                            +             &
               dvr_mat(5,10) * v(10,:)    &
                            +             &
               dvr_mat(5,11) * v(11,:)    

  v_scr(6,:) = v_scr(6,:)                 &
                            +             &  
               dvr_mat(6,1) * v(1,:)      &
                            +             &
               dvr_mat(6,2) * v(2,:)      &
                            +             &
               dvr_mat(6,3) * v(3,:)      &
                            +             &   
               dvr_mat(6,4) * v(4,:)      &
                            +             &
               dvr_mat(6,5) * v(5,:)      &   
                            +             &
               dvr_mat(6,6) * v(6,:)      &   
                            +             &
               dvr_mat(6,7) * v(7,:)      &   
                            +             &
               dvr_mat(6,8) * v(8,:)      &   
                            +             &
               dvr_mat(6,9) * v(9,:)      &  
                            +             &
               dvr_mat(6,10) * v(10,:)    & 
                            +             &
               dvr_mat(6,11) * v(11,:)    

  v_scr(7,:) = v_scr(7,:)                 &
                            +             &
               dvr_mat(7,1) * v(1,:)      &
                            +             &
               dvr_mat(7,2) * v(2,:)      &
                            +             &
               dvr_mat(7,3) * v(3,:)      &
                            +             &   
               dvr_mat(7,4) * v(4,:)      &
                            +             &
               dvr_mat(7,5) * v(5,:)      &   
                            +             &
               dvr_mat(7,6) * v(6,:)      &   
                            +             &
               dvr_mat(7,7) * v(7,:)      &   
                            +             &
               dvr_mat(7,8) * v(8,:)      &   
                            +             &
               dvr_mat(7,9) * v(9,:)      &   
                            +             &
               dvr_mat(7,10) * v(10,:)    &
                            +             &
               dvr_mat(7,11) * v(11,:)    

  v_scr(8,:) = v_scr(8,:)                 &
                            +             &
               dvr_mat(8,1) * v(1,:)      &
                            +             &
               dvr_mat(8,2) * v(2,:)      &
                            +             &
               dvr_mat(8,3) * v(3,:)      &
                            +             &   
               dvr_mat(8,4) * v(4,:)      &
                            +             &
               dvr_mat(8,5) * v(5,:)      &  
                            +             &
               dvr_mat(8,6) * v(6,:)      &   
                            +             &
               dvr_mat(8,7) * v(7,:)      &   
                            +             &
               dvr_mat(8,8) * v(8,:)      &   
                            +             &
               dvr_mat(8,9) * v(9,:)      &   
                            +             &
               dvr_mat(8,10) * v(10,:)    &
                            +             &
               dvr_mat(8,11) * v(11,:)    

  v_scr(9,:) = v_scr(9,:)                 &
                            +             &    
               dvr_mat(9,1) * v(1,:)      &
                            +             &
               dvr_mat(9,2) * v(2,:)      &
                            +             &
               dvr_mat(9,3) * v(3,:)      &
                            +             &   
               dvr_mat(9,4) * v(4,:)      &
                            +             &
               dvr_mat(9,5) * v(5,:)      &   
                            +             &
               dvr_mat(9,6) * v(6,:)      &   
                            +             &
               dvr_mat(9,7) * v(7,:)      &   
                            +             &
               dvr_mat(9,8) * v(8,:)      &   
                            +             &
               dvr_mat(9,9) * v(9,:)      &   
                            +             &
               dvr_mat(9,10) * v(10,:)    &
                            +             &
               dvr_mat(9,11) * v(11,:)    

  v_scr(10,:) = v_scr(10,:)              &
                            +            &
               dvr_mat(10,1) * v(1,:)    &
                              +          &
               dvr_mat(10,2) * v(2,:)    &
                              +          &
               dvr_mat(10,3) * v(3,:)    &
                              +          &   
               dvr_mat(10,4) * v(4,:)    &
                              +          &
               dvr_mat(10,5) * v(5,:)    &   
                              +          &
               dvr_mat(10,6) * v(6,:)    &   
                              +          &
               dvr_mat(10,7) * v(7,:)    &   
                              +          &
               dvr_mat(10,8) * v(8,:)    &   
                              +          &
               dvr_mat(10,9) * v(9,:)    &   
                              +          &
               dvr_mat(10,10) * v(10,:)   & 
                            +             &
               dvr_mat(10,11) * v(11,:)    

  v_scr(11,:) = v_scr(11,:)              &
                            +            &
               dvr_mat(11,1) * v(1,:)    &
                              +          &
               dvr_mat(11,2) * v(2,:)    &
                              +          &
               dvr_mat(11,3) * v(3,:)    &
                              +          &   
               dvr_mat(11,4) * v(4,:)    &
                              +          &
               dvr_mat(11,5) * v(5,:)    &   
                              +          &
               dvr_mat(11,6) * v(6,:)    &   
                              +          &
               dvr_mat(11,7) * v(7,:)    &   
                              +          &
               dvr_mat(11,8) * v(8,:)    &   
                              +          &
               dvr_mat(11,9) * v(9,:)    &   
                              +          &
               dvr_mat(11,10) * v(10,:)  & 
                            +             &
               dvr_mat(11,11) * v(11,:)    
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_11_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_gen_z
!***begin prologue     v_mat_v_gen_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propator matrix vector multiplies
!***                   for a general FEDVR Hamiltonian.
!***                   The parameter nj is a dummy and can take on values
!***                   consistent with any dimensional problem.  In a 1D
!***                   case nj=1, in 2D nj=nx and in 3D nj=ny*nx.
!
!***references
!***routines called    ebcxx, ambcxx, apbcxx
!***end prologue       
!
  SUBROUTINE v_mat_v_gen_z(v,             &
                           v_scr,         &
                           dvr_mat,       &
                           ni,nj,nk)
  IMPLICIT NONE
  INTEGER                              :: ni, nj, nk
  INTEGER                              :: i, k
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(nk,nk)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   DO i=1,nk
      DO k=1,nk
         v_scr(k,:) = v_scr(k,:) + dvr_mat(k,i) * v(i,:) 
      END DO
   END DO
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE v_mat_v_gen_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_2_z
!***begin prologue     v_mat_v_2_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies for
!***                   the special case of 2*2 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_2_z
!
  SUBROUTINE v_mat_v_2_z(v,         &
                         v_scr,     &
                         dvr_mat,   &
                         ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(2,2)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = v_scr(1,:)              &
                            +          &
               dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   

  v_scr(2,:) = v_scr(2,:)              &
                            +          &
               dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_2_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_3_z
!***begin prologue     v_mat_v_3_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_3_z
!
  SUBROUTINE v_mat_v_3_z(v,            &
                         v_scr,        &
                         dvr_mat,      &
                         ni,nj)
  IMPLICIT NONE
  INTEGER                                :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)           :: v
  COMPLEX*16, DIMENSION(ni,nj)           :: v_scr
  REAL*8, DIMENSION(3,3)                 :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = v_scr(1,:)              &
                            +          &
               dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   

  v_scr(2,:) = v_scr(2,:)              &
                            +          &
               dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   

  v_scr(3,:) = v_scr(3,:)              &
                            +          &
               dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_3_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_4_z
!***begin prologue     v_mat_v_4_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 4*4 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_4_z
!
  SUBROUTINE v_mat_v_4_z(v,             &
                         v_scr,         &
                         dvr_mat,       &
                         ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(4,4)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = v_scr(1,:)              &
                            +          &
               dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   &
                            +          &   
               dvr_mat(1,4) * v(4,:)   

  v_scr(2,:) = v_scr(2,:)              &
                            +          &
               dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   &
                            +          &
               dvr_mat(2,4) * v(4,:)   

  v_scr(3,:) = v_scr(3,:)              &
                            +          &
               dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   &
                            +          &
               dvr_mat(3,4) * v(4,:)   

  v_scr(4,:) = v_scr(4,:)              &
                            +          &
               dvr_mat(4,1) * v(1,:)   &
                            +          &
               dvr_mat(4,2) * v(2,:)   &
                            +          &
               dvr_mat(4,3) * v(3,:)   &
                            +          &
               dvr_mat(4,4) * v(4,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_4_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_5_z
!***begin prologue     v_mat_v_5_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_5_z
!
  SUBROUTINE v_mat_v_5_z(v,          &
                         v_scr,      &
                         dvr_mat,    &
                         ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(5,5)               :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = v_scr(1,:)              &
                            +          &
               dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   &
                            +          &   
               dvr_mat(1,4) * v(4,:)   &
                            +          &
               dvr_mat(1,5) * v(5,:)    

  v_scr(2,:) = v_scr(2,:)              &
                            +          &
               dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   &
                            +          &
               dvr_mat(2,4) * v(4,:)   &
                            +          &
               dvr_mat(2,5) * v(5,:)   

  v_scr(3,:) = v_scr(3,:)              &
                            +          &
               dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   &
                            +          &
               dvr_mat(3,4) * v(4,:)   &
                            +          &
               dvr_mat(3,5) * v(5,:)   
  v_scr(4,:) = v_scr(4,:)              &
                            +          &
               dvr_mat(4,1) * v(1,:)   &
                            +          &
               dvr_mat(4,2) * v(2,:)   &
                            +          &
               dvr_mat(4,3) * v(3,:)   &
                            +          &
               dvr_mat(4,4) * v(4,:)   &
                            +          &
               dvr_mat(4,5) * v(5,:)    
  v_scr(5,:) = v_scr(5,:)              &
                            +          &
               dvr_mat(5,1) * v(1,:)   &
                            +          &
               dvr_mat(5,2) * v(2,:)   &
                            +          &
               dvr_mat(5,3) * v(3,:)   &
                            +          &
               dvr_mat(5,4) * v(4,:)   &
                            +          &
               dvr_mat(5,5) * v(5,:)    
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_5_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_6_z
!***begin prologue     v_mat_v_6_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_6_z
!
  SUBROUTINE v_mat_v_6_z(v,           &    
                         v_scr,       &
                         dvr_mat,     &
                         ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(6,6)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = v_scr(1,:)              &
                            +          &
               dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   &
                            +          &   
               dvr_mat(1,4) * v(4,:)   &
                            +          &
               dvr_mat(1,5) * v(5,:)   &   
                            +          &
               dvr_mat(1,6) * v(6,:)    

  v_scr(2,:) = v_scr(2,:)              &
                            +          &   
               dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   &
                            +          &   
               dvr_mat(2,4) * v(4,:)   &
                            +          &
               dvr_mat(2,5) * v(5,:)   &   
                            +          &
               dvr_mat(2,6) * v(6,:)      

  v_scr(3,:) = v_scr(3,:)              &
                            +          & 
               dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   &
                            +          &   
               dvr_mat(3,4) * v(4,:)   &
                            +          &
               dvr_mat(3,5) * v(5,:)   &   
                            +          &
               dvr_mat(3,6) * v(6,:)      

  v_scr(4,:) = v_scr(4,:)              &
                            +          &  
               dvr_mat(4,1) * v(1,:)   &
                            +          &
               dvr_mat(4,2) * v(2,:)   &
                            +          &
               dvr_mat(4,3) * v(3,:)   &
                            +          &   
               dvr_mat(4,4) * v(4,:)   &
                            +          &
               dvr_mat(4,5) * v(5,:)   &   
                            +          &
               dvr_mat(4,6) * v(6,:)      

  v_scr(5,:) = v_scr(5,:)              &
                            +          &  
               dvr_mat(5,1) * v(1,:)   &
                            +          &
               dvr_mat(5,2) * v(2,:)   &
                            +          &
               dvr_mat(5,3) * v(3,:)   &
                            +          &   
               dvr_mat(5,4) * v(4,:)   &
                            +          &
               dvr_mat(5,5) * v(5,:)   &   
                            +          &
               dvr_mat(5,6) * v(6,:)      

  v_scr(6,:) = v_scr(6,:)              &
                            +          &  
               dvr_mat(6,1) * v(1,:)   &
                            +          &
               dvr_mat(6,2) * v(2,:)   &
                            +          &
               dvr_mat(6,3) * v(3,:)   &
                            +          &   
               dvr_mat(6,4) * v(4,:)   &
                            +          &
               dvr_mat(6,5) * v(5,:)   &   
                            +          &
               dvr_mat(6,6) * v(6,:)      
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_6_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_7_z
!***begin prologue     v_mat_v_7_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_7_z
!
  SUBROUTINE v_mat_v_7_z(v,           &    
                         v_scr,       &
                         dvr_mat,     &
                         ni,nj)
                     
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(7,7)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = v_scr(1,:)              &
                            +          &
               dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   &
                            +          &   
               dvr_mat(1,4) * v(4,:)   &
                            +          &
               dvr_mat(1,5) * v(5,:)   &   
                            +          &
               dvr_mat(1,6) * v(6,:)   &   
                            +          &
               dvr_mat(1,7) * v(7,:)      

  v_scr(2,:) = v_scr(2,:)              &
                            +          &
               dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   &
                            +          &   
               dvr_mat(2,4) * v(4,:)   &
                            +          &
               dvr_mat(2,5) * v(5,:)   &   
                            +          &
               dvr_mat(2,6) * v(6,:)   &   
                            +          &
               dvr_mat(2,7) * v(7,:)    

  v_scr(3,:) = v_scr(3,:)              &
                            +          &  
               dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   &
                            +          &   
               dvr_mat(3,4) * v(4,:)   &
                            +          &
               dvr_mat(3,5) * v(5,:)   &   
                            +          &
               dvr_mat(3,6) * v(6,:)   &   
                            +          &
               dvr_mat(3,7) * v(7,:)    

  v_scr(4,:) = v_scr(4,:)              &
                            +          &  
               dvr_mat(4,1) * v(1,:)   &
                            +          &
               dvr_mat(4,2) * v(2,:)   &
                            +          &
               dvr_mat(4,3) * v(3,:)   &
                            +          &   
               dvr_mat(4,4) * v(4,:)   &
                            +          &
               dvr_mat(4,5) * v(5,:)   &   
                            +          &
               dvr_mat(4,6) * v(6,:)   &   
                            +          &
                dvr_mat(4,7) * v(7,:)   

  v_scr(5,:) = v_scr(5,:)              &
                            +          &   
               dvr_mat(5,1) * v(1,:)   &
                            +          &
               dvr_mat(5,2) * v(2,:)   &
                            +          &
               dvr_mat(5,3) * v(3,:)   &
                            +          &   
               dvr_mat(5,4) * v(4,:)   &
                            +          &
               dvr_mat(5,5) * v(5,:)   &   
                            +          &
               dvr_mat(5,6) * v(6,:)   &   
                            +          &
               dvr_mat(5,7) * v(7,:)    

  v_scr(6,:) = v_scr(6,:)              &
                            +          &  
               dvr_mat(6,1) * v(1,:)   &
                            +          &
               dvr_mat(6,2) * v(2,:)   &
                            +          &
               dvr_mat(6,3) * v(3,:)   &
                            +          &   
               dvr_mat(6,4) * v(4,:)   &
                            +          &
               dvr_mat(6,5) * v(5,:)   &   
                            +          &
               dvr_mat(6,6) * v(6,:)   &   
                            +          &
               dvr_mat(6,7) * v(7,:)    

  v_scr(7,:) = v_scr(7,:)              &
                            +          &  
               dvr_mat(7,1) * v(1,:)   &
                            +          &
               dvr_mat(7,2) * v(2,:)   &
                            +          &
               dvr_mat(7,3) * v(3,:)   &
                            +          &   
               dvr_mat(7,4) * v(4,:)   &
                            +          &
               dvr_mat(7,5) * v(5,:)   &   
                            +          &
               dvr_mat(7,6) * v(6,:)   &   
                            +          &
               dvr_mat(7,7) * v(7,:)      
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_7_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_8_z
!***begin prologue     v_mat_v_8_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_8_z
!
  SUBROUTINE v_mat_v_8_z(v,             &
                         v_scr,         &
                         dvr_mat,       &
                         ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(8,8)               :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  v_scr(1,:) = v_scr(1,:)                 &
                            +             &
               dvr_mat(1,1) * v(1,:)      &           
                            +             &
               dvr_mat(1,2) * v(2,:)      &
                            +             &
               dvr_mat(1,3) * v(3,:)      &
                            +             &   
               dvr_mat(1,4) * v(4,:)      &
                            +             &
               dvr_mat(1,5) * v(5,:)      &   
                            +             &
               dvr_mat(1,6) * v(6,:)      &   
                            +             &
               dvr_mat(1,7) * v(7,:)      &   
                            +             &
               dvr_mat(1,8) * v(8,:)      

  v_scr(2,:) = v_scr(2,:)                 &
                               +          &
               dvr_mat(2,1) * v(1,:)      &
                            +             &
               dvr_mat(2,2) * v(2,:)      &
                            +             &
               dvr_mat(2,3) * v(3,:)      &
                            +             &   
               dvr_mat(2,4) * v(4,:)      &
                            +             &
               dvr_mat(2,5) * v(5,:)      &   
                            +             &
               dvr_mat(2,6) * v(6,:)      &   
                            +             &
               dvr_mat(2,7) * v(7,:)      &   
                            +             &
               dvr_mat(2,8) * v(8,:) 
  v_scr(3,:) = v_scr(3,:)                 &
                            +             &     
               dvr_mat(3,1) * v(1,:)      &
                            +             &
               dvr_mat(3,2) * v(2,:)      &
                            +             &
               dvr_mat(3,3) * v(3,:)      &
                            +             &   
               dvr_mat(3,4) * v(4,:)      &
                            +             &
               dvr_mat(3,5) * v(5,:)      &   
                            +             &
               dvr_mat(3,6) * v(6,:)      &   
                            +             &
               dvr_mat(3,7) * v(7,:)      &   
                            +             &
               dvr_mat(3,8) * v(8,:)      

  v_scr(4,:) = v_scr(4,:)                 &
                            +             &
               dvr_mat(4,1) * v(1,:)      &
                            +             &
               dvr_mat(4,2) * v(2,:)      &
                            +             &
               dvr_mat(4,3) * v(3,:)      &
                            +             &   
               dvr_mat(4,4) * v(4,:)      &
                            +             &
               dvr_mat(4,5) * v(5,:)      &   
                            +             &
               dvr_mat(4,6) * v(6,:)      &   
                            +             &
               dvr_mat(4,7) * v(7,:)      &   
                            +             &
               dvr_mat(4,8) * v(8,:)      

  v_scr(5,:) = v_scr(5,:)                &
                            +             &
               dvr_mat(5,1) * v(1,:)      &
                            +             &
               dvr_mat(5,2) * v(2,:)      &
                            +             &
               dvr_mat(5,3) * v(3,:)      &
                            +             &   
               dvr_mat(5,4) * v(4,:)      &
                            +             &
               dvr_mat(5,5) * v(5,:)      &   
                            +             &
               dvr_mat(5,6) * v(6,:)      &   
                            +             &
               dvr_mat(5,7) * v(7,:)      &   
                            +             &
               dvr_mat(5,8) * v(8,:)     
   
  v_scr(6,:) = v_scr(6,:)                 &
                            +             &
               dvr_mat(6,1) * v(1,:)      &
                            +             &
               dvr_mat(6,2) * v(2,:)      &
                            +             &
               dvr_mat(6,3) * v(3,:)      &
                            +             &   
               dvr_mat(6,4) * v(4,:)      &
                            +             &
               dvr_mat(6,5) * v(5,:)      &   
                            +             &
               dvr_mat(6,6) * v(6,:)      &   
                            +             &
               dvr_mat(6,7) * v(7,:)      &   
                            +             &
               dvr_mat(6,8) * v(8,:)        

  v_scr(7,:) = v_scr(7,:)                 &
                            +             &
               dvr_mat(7,1) * v(1,:)      &
                            +             &
               dvr_mat(7,2) * v(2,:)      &
                            +             &
               dvr_mat(7,3) * v(3,:)      &
                            +             &   
               dvr_mat(7,4) * v(4,:)      &
                            +             &
               dvr_mat(7,5) * v(5,:)      &   
                            +             &
               dvr_mat(7,6) * v(6,:)      &   
                            +             &
               dvr_mat(7,7) * v(7,:)      &   
                            +             &
               dvr_mat(7,8) * v(8,:)        

  v_scr(8,:) = v_scr(8,:)                 &
                            +             &
               dvr_mat(8,1) * v(1,:)      &
                            +             &
               dvr_mat(8,2) * v(2,:)      &
                            +             &
               dvr_mat(8,3) * v(3,:)      &
                            +             &   
               dvr_mat(8,4) * v(4,:)      &
                            +             &
               dvr_mat(8,5) * v(5,:)      &  
                            +             &
               dvr_mat(8,6) * v(6,:)      &   
                            +             &
               dvr_mat(8,7) * v(7,:)      &   
                            +             &
               dvr_mat(8,8) * v(8,:)        
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_8_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_9_z
!***begin prologue     v_mat_v_9_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_9_z
!
  SUBROUTINE v_mat_v_9_z(v,           &
                         v_scr,       &
                         dvr_mat,     &
                         ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16 , DIMENSION(ni,nj)        :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(9,9)               :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = v_scr(1,:)                 &
                            +             &
               dvr_mat(1,1) * v(1,:)      &
                            +             &
               dvr_mat(1,2) * v(2,:)      &
                            +             &
               dvr_mat(1,3) * v(3,:)      &
                            +             &   
               dvr_mat(1,4) * v(4,:)      &
                            +             &
               dvr_mat(1,5) * v(5,:)      &   
                            +             &
               dvr_mat(1,6) * v(6,:)      &   
                            +             &
               dvr_mat(1,7) * v(7,:)      &   
                            +             &
               dvr_mat(1,8) * v(8,:)      &   
                            +             &
               dvr_mat(1,9) * v(9,:)    

  v_scr(2,:) = v_scr(2,:)                 &
                            +             &  
               dvr_mat(2,1) * v(1,:)      &
                            +             &
               dvr_mat(2,2) * v(2,:)      &
                            +             &
               dvr_mat(2,3) * v(3,:)      &
                            +             &   
               dvr_mat(2,4) * v(4,:)      &
                            +             &
               dvr_mat(2,5) * v(5,:)      &
                            +             &
               dvr_mat(2,6) * v(6,:)      &   
                            +             &
               dvr_mat(2,7) * v(7,:)      &   
                            +             &
               dvr_mat(2,8) * v(8,:)      &   
                            +             &
               dvr_mat(2,9) * v(9,:)    

  v_scr(3,:) = v_scr(3,:)                 &
                            +             &  
               dvr_mat(3,1) * v(1,:)      &
                            +             &
               dvr_mat(3,2) * v(2,:)      &
                            +             &
               dvr_mat(3,3) * v(3,:)      &
                            +             &   
               dvr_mat(3,4) * v(4,:)      &
                            +             &
               dvr_mat(3,5) * v(5,:)      &   
                            +             &
               dvr_mat(3,6) * v(6,:)      &   
                            +             &
               dvr_mat(3,7) * v(7,:)      &   
                            +             &
               dvr_mat(3,8) * v(8,:)      &   
                            +             &
               dvr_mat(3,9) * v(9,:)      


  v_scr(4,:) = v_scr(4,:)                 &
                            +             &
               dvr_mat(4,1) * v(1,:)      &
                            +             &
               dvr_mat(4,2) * v(2,:)      &
                            +             &
               dvr_mat(4,3) * v(3,:)      &
                            +             &   
               dvr_mat(4,4) * v(4,:)      &
                            +             &
               dvr_mat(4,5) * v(5,:)      &   
                            +             &
               dvr_mat(4,6) * v(6,:)      &   
                            +             &
               dvr_mat(4,7) * v(7,:)      &   
                            +             &
               dvr_mat(4,8) * v(8,:)      &   
                            +             &
               dvr_mat(4,9) * v(9,:)      

  v_scr(5,:) = v_scr(5,:)                 &
                            +             &
               dvr_mat(5,1) * v(1,:)      &
                            +             &
               dvr_mat(5,2) * v(2,:)      &
                            +             &
               dvr_mat(5,3) * v(3,:)      &
                            +             &   
               dvr_mat(5,4) * v(4,:)      &
                            +             &
               dvr_mat(5,5) * v(5,:)      &   
                            +             &
               dvr_mat(5,6) * v(6,:)      &   
                            +             &
               dvr_mat(5,7) * v(7,:)      &   
                            +             &
               dvr_mat(5,8) * v(8,:)      &   
                            +             &
               dvr_mat(5,9) * v(9,:)      

  v_scr(6,:) = v_scr(6,:)                 &
                            +             &
               dvr_mat(6,1) * v(1,:)      &
                            +             &
               dvr_mat(6,2) * v(2,:)      &
                            +             &
               dvr_mat(6,3) * v(3,:)      &
                            +             &   
               dvr_mat(6,4) * v(4,:)      &
                            +             &
               dvr_mat(6,5) * v(5,:)      &   
                            +             &
               dvr_mat(6,6) * v(6,:)      &   
                            +             &
               dvr_mat(6,7) * v(7,:)      &   
                            +             &
               dvr_mat(6,8) * v(8,:)      &   
                            +             &
               dvr_mat(6,9) * v(9,:)      

  v_scr(7,:) = v_scr(7,:)                 &
                            +             &
               dvr_mat(7,1) * v(1,:)      &
                            +             &
               dvr_mat(7,2) * v(2,:)      &
                            +             &
               dvr_mat(7,3) * v(3,:)      &
                            +             &   
               dvr_mat(7,4) * v(4,:)      &
                            +             &
               dvr_mat(7,5) * v(5,:)      &   
                            +             &
               dvr_mat(7,6) * v(6,:)      &  
                            +             &
               dvr_mat(7,7) * v(7,:)      &   
                            +             &
               dvr_mat(7,8) * v(8,:)      &   
                            +             &
               dvr_mat(7,9) * v(9,:)      

  v_scr(8,:) = v_scr(8,:)                 &
                            +             &
               dvr_mat(8,1) * v(1,:)      &
                            +             &
               dvr_mat(8,2) * v(2,:)      &
                            +             &
               dvr_mat(8,3) * v(3,:)      &
                            +             &   
               dvr_mat(8,4) * v(4,:)      &
                            +             &
               dvr_mat(8,5) * v(5,:)      &   
                            +             &
               dvr_mat(8,6) * v(6,:)      &   
                            +             &
               dvr_mat(8,7) * v(7,:)      &   
                            +             &
               dvr_mat(8,8) * v(8,:)      &   
                            +             &
               dvr_mat(8,9) * v(9,:)      

  v_scr(9,:) = v_scr(9,:)                 &
                            +             &
               dvr_mat(9,1) * v(1,:)      &
                            +             &
               dvr_mat(9,2) * v(2,:)      &
                            +             &
               dvr_mat(9,3) * v(3,:)      &
                            +             &   
               dvr_mat(9,4) * v(4,:)      &
                            +             &
               dvr_mat(9,5) * v(5,:)      &   
                            +             &
               dvr_mat(9,6) * v(6,:)      &   
                            +             &
               dvr_mat(9,7) * v(7,:)      &   
                            +             &
               dvr_mat(9,8) * v(8,:)      &   
                            +             &
               dvr_mat(9,9) * v(9,:)      
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_9_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_10_z
!***begin prologue     v_mat_v_10_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_10_z
!
  SUBROUTINE v_mat_v_10_z(v,          &
                          v_scr,      &
                          dvr_mat,    &
                          ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(10,10)             :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = v_scr(1,:)                 &
                            +             &
               dvr_mat(1,1) * v(1,:)      &
                            +             &
               dvr_mat(1,2) * v(2,:)      &
                            +             &
               dvr_mat(1,3) * v(3,:)      &
                            +             &   
               dvr_mat(1,4) * v(4,:)      &
                            +             &
               dvr_mat(1,5) * v(5,:)      &   
                            +             &
               dvr_mat(1,6) * v(6,:)      &   
                            +             &
               dvr_mat(1,7) * v(7,:)      &   
                            +             &
               dvr_mat(1,8) * v(8,:)      &   
                            +             &
               dvr_mat(1,9) * v(9,:)      &   
                            +             &
               dvr_mat(1,10) * v(10,:)    

  v_scr(2,:) = v_scr(2,:)                 &
                            +             &
               dvr_mat(2,1) * v(1,:)      &
                            +             &
               dvr_mat(2,2) * v(2,:)      &
                            +             &
               dvr_mat(2,3) * v(3,:)      &
                            +             &   
               dvr_mat(2,4) * v(4,:)      &
                            +             &
               dvr_mat(2,5) * v(5,:)      &   
                            +             &
               dvr_mat(2,6) * v(6,:)      &   
                            +             &
               dvr_mat(2,7) * v(7,:)      &   
                            +             &
               dvr_mat(2,8) * v(8,:)      &   
                            +             &
               dvr_mat(2,9) * v(9,:)      &   
                            +             &
               dvr_mat(2,10) * v(10,:)    

  v_scr(3,:) = v_scr(3,:)                 &
                            +             &
               dvr_mat(3,1) * v(1,:)      &
                            +             &
               dvr_mat(3,2) * v(2,:)      &
                            +             &
               dvr_mat(3,3) * v(3,:)      &
                            +             &   
               dvr_mat(3,4) * v(4,:)      &
                            +             &
               dvr_mat(3,5) * v(5,:)      &   
                            +             &
               dvr_mat(3,6) * v(6,:)      &   
                            +             &
               dvr_mat(3,7) * v(7,:)      &   
                            +             &
               dvr_mat(3,8) * v(8,:)      &   
                            +             &
               dvr_mat(3,9) * v(9,:)      &   
                            +             &
               dvr_mat(3,10) * v(10,:)    

  v_scr(4,:) = v_scr(4,:)                 &
                            +             &
               dvr_mat(4,1) * v(1,:)      &
                            +             &
               dvr_mat(4,2) * v(2,:)      &
                            +             &
               dvr_mat(4,3) * v(3,:)      &
                            +             &   
               dvr_mat(4,4) * v(4,:)      &
                            +             &
               dvr_mat(4,5) * v(5,:)      &   
                            +             &
               dvr_mat(4,6) * v(6,:)      &   
                            +             &
               dvr_mat(4,7) * v(7,:)      &   
                            +             &
               dvr_mat(4,8) * v(8,:)      &   
                            +             &
               dvr_mat(4,9) * v(9,:)      &   
                            +             &
               dvr_mat(4,10) * v(10,:)    

  v_scr(5,:) = v_scr(5,:)                 &
                            +             &
               dvr_mat(5,1) * v(1,:)      &
                            +             &
               dvr_mat(5,2) * v(2,:)      &
                            +             &
               dvr_mat(5,3) * v(3,:)      &
                            +             &   
               dvr_mat(5,4) * v(4,:)      &
                            +             &
               dvr_mat(5,5) * v(5,:)      &   
                            +             &
               dvr_mat(5,6) * v(6,:)      &   
                            +             &
               dvr_mat(5,7) * v(7,:)      &   
                            +             &
               dvr_mat(5,8) * v(8,:)      &   
                            +             &
               dvr_mat(5,9) * v(9,:)      &   
                            +             &
               dvr_mat(5,10) * v(10,:)  

  v_scr(6,:) = v_scr(6,:)                 &
                            +             &  
               dvr_mat(6,1) * v(1,:)      &
                            +             &
               dvr_mat(6,2) * v(2,:)      &
                            +             &
               dvr_mat(6,3) * v(3,:)      &
                            +             &   
               dvr_mat(6,4) * v(4,:)      &
                            +             &
               dvr_mat(6,5) * v(5,:)      &   
                            +             &
               dvr_mat(6,6) * v(6,:)      &   
                            +             &
               dvr_mat(6,7) * v(7,:)      &   
                            +             &
               dvr_mat(6,8) * v(8,:)      &   
                            +             &
               dvr_mat(6,9) * v(9,:)      &  
                            +             &
               dvr_mat(6,10) * v(10,:)     

  v_scr(7,:) = v_scr(7,:)                 &
                            +             &
               dvr_mat(7,1) * v(1,:)      &
                            +             &
               dvr_mat(7,2) * v(2,:)      &
                            +             &
               dvr_mat(7,3) * v(3,:)      &
                            +             &   
               dvr_mat(7,4) * v(4,:)      &
                            +             &
               dvr_mat(7,5) * v(5,:)      &   
                            +             &
               dvr_mat(7,6) * v(6,:)      &   
                            +             &
               dvr_mat(7,7) * v(7,:)      &   
                            +             &
               dvr_mat(7,8) * v(8,:)      &   
                            +             &
               dvr_mat(7,9) * v(9,:)      &   
                            +             &
               dvr_mat(7,10) * v(10,:)   

  v_scr(8,:) = v_scr(8,:)                 &
                            +             &
               dvr_mat(8,1) * v(1,:)      &
                            +             &
               dvr_mat(8,2) * v(2,:)      &
                            +             &
               dvr_mat(8,3) * v(3,:)      &
                            +             &   
               dvr_mat(8,4) * v(4,:)      &
                            +             &
               dvr_mat(8,5) * v(5,:)      &  
                            +             &
               dvr_mat(8,6) * v(6,:)      &   
                            +             &
               dvr_mat(8,7) * v(7,:)      &   
                            +             &
               dvr_mat(8,8) * v(8,:)      &   
                            +             &
               dvr_mat(8,9) * v(9,:)      &   
                            +             &
               dvr_mat(8,10) * v(10,:)  

  v_scr(9,:) = v_scr(9,:)                 &
                            +             &    
               dvr_mat(9,1) * v(1,:)      &
                            +             &
               dvr_mat(9,2) * v(2,:)      &
                            +             &
               dvr_mat(9,3) * v(3,:)      &
                            +             &   
               dvr_mat(9,4) * v(4,:)      &
                            +             &
               dvr_mat(9,5) * v(5,:)      &   
                            +             &
               dvr_mat(9,6) * v(6,:)      &   
                            +             &
               dvr_mat(9,7) * v(7,:)      &   
                            +             &
               dvr_mat(9,8) * v(8,:)      &   
                            +             &
               dvr_mat(9,9) * v(9,:)      &   
                            +             &
               dvr_mat(9,10) * v(10,:)    

  v_scr(10,:) = v_scr(10,:)              &
                            +            &
               dvr_mat(10,1) * v(1,:)    &
                              +          &
               dvr_mat(10,2) * v(2,:)    &
                              +          &
               dvr_mat(10,3) * v(3,:)    &
                              +          &   
               dvr_mat(10,4) * v(4,:)    &
                              +          &
               dvr_mat(10,5) * v(5,:)    &   
                              +          &
               dvr_mat(10,6) * v(6,:)    &   
                              +          &
               dvr_mat(10,7) * v(7,:)    &   
                              +          &
               dvr_mat(10,8) * v(8,:)    &   
                              +          &
               dvr_mat(10,9) * v(9,:)    &   
                              +          &
               dvr_mat(10,10) * v(10,:)    
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_10_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_mat_v_11_z
!***begin prologue     v_mat_v_11_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_mat_v_11_z
!
  SUBROUTINE v_mat_v_11_z(v,          &
                          v_scr,      &
                          dvr_mat,    &
                          ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(11,11)             :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = v_scr(1,:)                 &
                            +             &
               dvr_mat(1,1) * v(1,:)      &
                            +             &
               dvr_mat(1,2) * v(2,:)      &
                            +             &
               dvr_mat(1,3) * v(3,:)      &
                            +             &   
               dvr_mat(1,4) * v(4,:)      &
                            +             &
               dvr_mat(1,5) * v(5,:)      &   
                            +             &
               dvr_mat(1,6) * v(6,:)      &   
                            +             &
               dvr_mat(1,7) * v(7,:)      &   
                            +             &
               dvr_mat(1,8) * v(8,:)      &   
                            +             &
               dvr_mat(1,9) * v(9,:)      &   
                            +             &
               dvr_mat(1,10) * v(10,:)    &
                            +             &
               dvr_mat(1,11) * v(11,:)    

  v_scr(2,:) = v_scr(2,:)                 &
                            +             &
               dvr_mat(2,1) * v(1,:)      &
                            +             &
               dvr_mat(2,2) * v(2,:)      &
                            +             &
               dvr_mat(2,3) * v(3,:)      &
                            +             &   
               dvr_mat(2,4) * v(4,:)      &
                            +             &
               dvr_mat(2,5) * v(5,:)      &   
                            +             &
               dvr_mat(2,6) * v(6,:)      &   
                            +             &
               dvr_mat(2,7) * v(7,:)      &   
                            +             &
               dvr_mat(2,8) * v(8,:)      &   
                            +             &
               dvr_mat(2,9) * v(9,:)      &   
                            +             &
               dvr_mat(2,10) * v(10,:)    &
                            +             &
               dvr_mat(2,11) * v(11,:)    

  v_scr(3,:) = v_scr(3,:)                 &
                            +             &
               dvr_mat(3,1) * v(1,:)      &
                            +             &
               dvr_mat(3,2) * v(2,:)      &
                            +             &
               dvr_mat(3,3) * v(3,:)      &
                            +             &   
               dvr_mat(3,4) * v(4,:)      &
                            +             &
               dvr_mat(3,5) * v(5,:)      &   
                            +             &
               dvr_mat(3,6) * v(6,:)      &   
                            +             &
               dvr_mat(3,7) * v(7,:)      &   
                            +             &
               dvr_mat(3,8) * v(8,:)      &   
                            +             &
               dvr_mat(3,9) * v(9,:)      &   
                            +             &
               dvr_mat(3,10) * v(10,:)    &
                            +             &
               dvr_mat(3,11) * v(11,:)    

  v_scr(4,:) = v_scr(4,:)                 &
                            +             &
               dvr_mat(4,1) * v(1,:)      &
                            +             &
               dvr_mat(4,2) * v(2,:)      &
                            +             &
               dvr_mat(4,3) * v(3,:)      &
                            +             &   
               dvr_mat(4,4) * v(4,:)      &
                            +             &
               dvr_mat(4,5) * v(5,:)      &   
                            +             &
               dvr_mat(4,6) * v(6,:)      &   
                            +             &
               dvr_mat(4,7) * v(7,:)      &   
                            +             &
               dvr_mat(4,8) * v(8,:)      &   
                            +             &
               dvr_mat(4,9) * v(9,:)      &   
                            +             &
               dvr_mat(4,10) * v(10,:)    &
                            +             &
               dvr_mat(4,11) * v(11,:)    

  v_scr(5,:) = v_scr(5,:)                 &
                            +             &
               dvr_mat(5,1) * v(1,:)      &
                            +             &
               dvr_mat(5,2) * v(2,:)      &
                            +             &
               dvr_mat(5,3) * v(3,:)      &
                            +             &   
               dvr_mat(5,4) * v(4,:)      &
                            +             &
               dvr_mat(5,5) * v(5,:)      &   
                            +             &
               dvr_mat(5,6) * v(6,:)      &   
                            +             &
               dvr_mat(5,7) * v(7,:)      &   
                            +             &
               dvr_mat(5,8) * v(8,:)      &   
                            +             &
               dvr_mat(5,9) * v(9,:)      &   
                            +             &
               dvr_mat(5,10) * v(10,:)    &
                            +             &
               dvr_mat(5,11) * v(11,:)    

  v_scr(6,:) = v_scr(6,:)                 &
                            +             &  
               dvr_mat(6,1) * v(1,:)      &
                            +             &
               dvr_mat(6,2) * v(2,:)      &
                            +             &
               dvr_mat(6,3) * v(3,:)      &
                            +             &   
               dvr_mat(6,4) * v(4,:)      &
                            +             &
               dvr_mat(6,5) * v(5,:)      &   
                            +             &
               dvr_mat(6,6) * v(6,:)      &   
                            +             &
               dvr_mat(6,7) * v(7,:)      &   
                            +             &
               dvr_mat(6,8) * v(8,:)      &   
                            +             &
               dvr_mat(6,9) * v(9,:)      &  
                            +             &
               dvr_mat(6,10) * v(10,:)    & 
                            +             &
               dvr_mat(6,11) * v(11,:)    

  v_scr(7,:) = v_scr(7,:)                 &
                            +             &
               dvr_mat(7,1) * v(1,:)      &
                            +             &
               dvr_mat(7,2) * v(2,:)      &
                            +             &
               dvr_mat(7,3) * v(3,:)      &
                            +             &   
               dvr_mat(7,4) * v(4,:)      &
                            +             &
               dvr_mat(7,5) * v(5,:)      &   
                            +             &
               dvr_mat(7,6) * v(6,:)      &   
                            +             &
               dvr_mat(7,7) * v(7,:)      &   
                            +             &
               dvr_mat(7,8) * v(8,:)      &   
                            +             &
               dvr_mat(7,9) * v(9,:)      &   
                            +             &
               dvr_mat(7,10) * v(10,:)    &
                            +             &
               dvr_mat(7,11) * v(11,:)    

  v_scr(8,:) = v_scr(8,:)                 &
                            +             &
               dvr_mat(8,1) * v(1,:)      &
                            +             &
               dvr_mat(8,2) * v(2,:)      &
                            +             &
               dvr_mat(8,3) * v(3,:)      &
                            +             &   
               dvr_mat(8,4) * v(4,:)      &
                            +             &
               dvr_mat(8,5) * v(5,:)      &  
                            +             &
               dvr_mat(8,6) * v(6,:)      &   
                            +             &
               dvr_mat(8,7) * v(7,:)      &   
                            +             &
               dvr_mat(8,8) * v(8,:)      &   
                            +             &
               dvr_mat(8,9) * v(9,:)      &   
                            +             &
               dvr_mat(8,10) * v(10,:)    &
                            +             &
               dvr_mat(8,11) * v(11,:)    

  v_scr(9,:) = v_scr(9,:)                 &
                            +             &    
               dvr_mat(9,1) * v(1,:)      &
                            +             &
               dvr_mat(9,2) * v(2,:)      &
                            +             &
               dvr_mat(9,3) * v(3,:)      &
                            +             &   
               dvr_mat(9,4) * v(4,:)      &
                            +             &
               dvr_mat(9,5) * v(5,:)      &   
                            +             &
               dvr_mat(9,6) * v(6,:)      &   
                            +             &
               dvr_mat(9,7) * v(7,:)      &   
                            +             &
               dvr_mat(9,8) * v(8,:)      &   
                            +             &
               dvr_mat(9,9) * v(9,:)      &   
                            +             &
               dvr_mat(9,10) * v(10,:)    &
                            +             &
               dvr_mat(9,11) * v(11,:)    

  v_scr(10,:) = v_scr(10,:)              &
                            +            &
               dvr_mat(10,1) * v(1,:)    &
                              +          &
               dvr_mat(10,2) * v(2,:)    &
                              +          &
               dvr_mat(10,3) * v(3,:)    &
                              +          &   
               dvr_mat(10,4) * v(4,:)    &
                              +          &
               dvr_mat(10,5) * v(5,:)    &   
                              +          &
               dvr_mat(10,6) * v(6,:)    &   
                              +          &
               dvr_mat(10,7) * v(7,:)    &   
                              +          &
               dvr_mat(10,8) * v(8,:)    &   
                              +          &
               dvr_mat(10,9) * v(9,:)    &   
                              +          &
               dvr_mat(10,10) * v(10,:)   & 
                            +             &
               dvr_mat(10,11) * v(11,:)    

  v_scr(11,:) = v_scr(11,:)              &
                            +            &
               dvr_mat(11,1) * v(1,:)    &
                              +          &
               dvr_mat(11,2) * v(2,:)    &
                              +          &
               dvr_mat(11,3) * v(3,:)    &
                              +          &   
               dvr_mat(11,4) * v(4,:)    &
                              +          &
               dvr_mat(11,5) * v(5,:)    &   
                              +          &
               dvr_mat(11,6) * v(6,:)    &   
                              +          &
               dvr_mat(11,7) * v(7,:)    &   
                              +          &
               dvr_mat(11,8) * v(8,:)    &   
                              +          &
               dvr_mat(11,9) * v(9,:)    &   
                              +          &
               dvr_mat(11,10) * v(10,:)  & 
                            +             &
               dvr_mat(11,11) * v(11,:)    
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_mat_v_11_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**begin prologue     v_v_mat
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       v_v_mat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck v_v_mat_gen_d
!***begin prologue     v_v_mat_gen_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multipiles
!***                   for the general FEDVR HAmiltonian. 
!***                   The routine is not needed in 1D problems. 
!***                   In 2D, the i index runs over the number
!                      of points in the y coordinate while in 3D, it runs over
!                      the product of the number of points in the z and y coordinates.
!
!***references
!***routines called    ebcxx, ambcxx, apbcxx
!***end prologue   
!
  SUBROUTINE v_v_mat_gen_d(v,          &
                           v_scr,      &
                           dvr_mat,    & 
                           ni,nj,nk)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj, nk
  INTEGER                              :: j, k
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(nk,nk)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
!
!
  DO k=1,nk
     DO j=1,nk
        v_scr(:,j) = v_scr(:,j) + v(:,k) * dvr_mat(k,j) 
      END DO
   END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE v_v_mat_gen_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_2_d
!***begin prologue     v_v_mat_2_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagation matrix multiplies
!**                    for the special case of 2*2 matrices.
!
!***references
!***routines called
!***end prologue       v_v_mat_2_d
  SUBROUTINE v_v_mat_2_d(v,             &
                         v_scr,         &
                         dvr_mat,       & 
                         ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(2,2)               :: dvr_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  v_scr(:,1) = v_scr(:,1)                &
                       +                 &
               v(:,1)  * dvr_mat(1,1)    &
                       +                 &
               v(:,2)  * dvr_mat(2,1)    
 
  v_scr(:,2) = v_scr(:,2)                &
                       +                 &
               v(:,1)  * dvr_mat(1,2)    &
                       +                 &
               v(:,2)  * dvr_mat(2,2)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_2_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_3_d
!***begin prologue     v_v_mat_3_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_v_mat_3_d
!
  SUBROUTINE v_v_mat_3_d(v,          &
                         v_scr,      &
                         dvr_mat,    & 
                         ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(3,3)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       
  v_scr(:,1) = v_scr(:,1)                &
                       +                 &
               v(:,1)  * dvr_mat(1,1)    &
                       +                 &
               v(:,2)  * dvr_mat(2,1)    &
                       +                 &
               v(:,3)  * dvr_mat(3,1)    

  v_scr(:,2) = v_scr(:,2)                &
                       +                 &
               v(:,1)  * dvr_mat(1,2)    &
                       +                 &
               v(:,2)  * dvr_mat(2,2)    &
                       +                 &
               v(:,3)  * dvr_mat(3,2)    

  v_scr(:,3) = v_scr(:,3)                &
                       +                 &
               v(:,1)  * dvr_mat(1,3)    &
                       +                 &
               v(:,2)  * dvr_mat(2,3)    &
                       +                 &
               v(:,3)  * dvr_mat(3,3)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_3_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_4_d
!***begin prologue     v_v_mat_4_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_v_mat_4_d
!
  SUBROUTINE v_v_mat_4_d(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(4,4)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
  v_scr(:,1) = v_scr(:,1)                &
                       +                 &
               v(:,1)  * dvr_mat(1,1)    &
                       +                 &
               v(:,2)  * dvr_mat(2,1)    &
                       +                 &
               v(:,3)  * dvr_mat(3,1)    &
                       +                 &
               v(:,4)  * dvr_mat(4,1)    

  v_scr(:,2) = v_scr(:,2)                &
                       +                 &
               v(:,1)  * dvr_mat(1,2)    &
                       +                 &
               v(:,2)  * dvr_mat(2,2)    &
                       +                 &
               v(:,3)  * dvr_mat(3,2)    &
                       +                 &
               v(:,4)  * dvr_mat(4,2)    

  v_scr(:,3) = v_scr(:,3)                &
                       +                 &
               v(:,1)  * dvr_mat(1,3)    &
                       +                 &
               v(:,2)  * dvr_mat(2,3)    &
                       +                 &
               v(:,3)  * dvr_mat(3,3)    &
                       +                 &
               v(:,4)  * dvr_mat(4,3)    

  v_scr(:,4) = v_scr(:,4)                &
                       +                 &
               v(:,1)  * dvr_mat(1,4)    &
                       +                 &
               v(:,2)  * dvr_mat(2,4)    &
                       +                 &
               v(:,3)  * dvr_mat(3,4)    &
                       +                 &
               v(:,4)  * dvr_mat(4,4)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_4_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_5_d
!***begin prologue     v_v_mat_5_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_mat_5_d
!
  SUBROUTINE v_v_mat_5_d(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(5,5)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_5_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_6_d
!***begin prologue     v_v_mat_6_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_mat_6_d
!
  SUBROUTINE v_v_mat_6_d(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(6,6)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)    

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_6_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_7_d
!***begin prologue     v_v_mat_7_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_mat_7_d
!
  SUBROUTINE v_v_mat_7_d(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(7,7)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)    

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)    

  v_scr(:,7) = v_scr(:,7)              &
                       +               &
               v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_7_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_8_d
!***begin prologue     v_v_mat_8_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_mat_8_d
!
  SUBROUTINE v_v_mat_8_d(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(8,8)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)  &
                       +               &
               v(:,8)  * dvr_mat(8,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)  &
                       +               &
               v(:,8)  * dvr_mat(8,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)  &
                       +               &
               v(:,8)  * dvr_mat(8,5)    

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)  &
                       +               &
               v(:,8)  * dvr_mat(8,6)    

  v_scr(:,7) = v_scr(:,7)              &
                       +               &
               v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)  &
                       +               &
               v(:,8)  * dvr_mat(8,7)    

  v_scr(:,8) = v_scr(:,8)              &
                       +               &
               v(:,1)  * dvr_mat(1,8)  &
                       +               &
               v(:,2)  * dvr_mat(2,8)  &
                       +               &
               v(:,3)  * dvr_mat(3,8)  &
                       +               &
               v(:,4)  * dvr_mat(4,8)  &
                       +               &
               v(:,5)  * dvr_mat(5,8)  &
                       +               &
               v(:,6)  * dvr_mat(6,8)  &
                       +               &
               v(:,7)  * dvr_mat(7,8)  &
                       +               &
               v(:,8)  * dvr_mat(8,8)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_8_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_9_d
!***begin prologue     v_v_mat_9_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_mat_9_d
!
  SUBROUTINE v_v_mat_9_d(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(9,9)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)  &
                       +               &
               v(:,8)  * dvr_mat(8,1)  &
                       +               &
               v(:,9)  * dvr_mat(9,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)  &
                       +               &
               v(:,9)  * dvr_mat(9,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)  &
                       +               &
               v(:,8)  * dvr_mat(8,3)  &
                       +               &
               v(:,9)  * dvr_mat(9,3)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)  &
                       +               &
               v(:,8)  * dvr_mat(8,4)  &
                       +               &
               v(:,9)  * dvr_mat(9,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)  &
                       +               &
               v(:,8)  * dvr_mat(8,5)  &
                       +               &
               v(:,9)  * dvr_mat(9,5)    

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)  &
                       +               &
               v(:,8)  * dvr_mat(8,6)  &
                       +               &
               v(:,9)  * dvr_mat(9,6)    

  v_scr(:,7) = v_scr(:,7)              &
                       +               &
               v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)  &
                       +               &
               v(:,8)  * dvr_mat(8,7)  &
                       +               &
               v(:,9)  * dvr_mat(9,7)    

  v_scr(:,8) = v_scr(:,8)              &
                       +               &
               v(:,1)  * dvr_mat(1,8)  &
                       +               &
               v(:,2)  * dvr_mat(2,8)  &
                       +               &
               v(:,3)  * dvr_mat(3,8)  &
                       +               &
               v(:,4)  * dvr_mat(4,8)  &
                       +               &
               v(:,5)  * dvr_mat(5,8)  &
                       +               &
               v(:,6)  * dvr_mat(6,8)  &
                       +               &
               v(:,7)  * dvr_mat(7,8)  &
                       +               &
               v(:,8)  * dvr_mat(8,8)  &
                       +               &
               v(:,9)  * dvr_mat(9,8)    

  v_scr(:,9) = v_scr(:,9)              &
                       +               &
               v(:,1)  * dvr_mat(1,9)  &
                       +               &
               v(:,2)  * dvr_mat(2,9)  &
                       +               &
               v(:,3)  * dvr_mat(3,9)  &
                       +               &
               v(:,4)  * dvr_mat(4,9)  &
                       +               &
               v(:,5)  * dvr_mat(5,9)  &
                       +               &
               v(:,6)  * dvr_mat(6,9)  &
                       +               &
               v(:,7)  * dvr_mat(7,9)  &
                       +               &
               v(:,8)  * dvr_mat(8,9)  &
                       +               &
               v(:,9)  * dvr_mat(9,9)    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_9_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_10_d
!***begin prologue     v_v_mat_10_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_mat_10_d
!
  SUBROUTINE v_v_mat_10_d(v,        &
                          v_scr,    &
                          dvr_mat,  & 
                          ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(10,10)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)  &
                       +               &
               v(:,8)  * dvr_mat(8,1)  &
                       +               &
               v(:,9)  * dvr_mat(9,1)  &
                       +               &
               v(:,10) * dvr_mat(10,1)   

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)  &
                       +               &
               v(:,9)  * dvr_mat(9,2)  &
                       +               &
               v(:,10) * dvr_mat(10,2)   

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)  &
                       +               &
               v(:,8)  * dvr_mat(8,3)  &
                       +               &
               v(:,9)  * dvr_mat(9,3)  &
                       +               &
               v(:,10) * dvr_mat(10,3)   

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)  &
                       +               &
               v(:,8)  * dvr_mat(8,4)  &
                       +               &
               v(:,9)  * dvr_mat(9,4)  &
                       +               &
               v(:,10) * dvr_mat(10,4)   

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)  &
                       +               &
               v(:,8)  * dvr_mat(8,5)  &
                       +               &
               v(:,9)  * dvr_mat(9,5)  &
                       +               &
               v(:,10) * dvr_mat(10,5)   

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)  &
                       +               &
               v(:,8)  * dvr_mat(8,6)  &
                       +               &
               v(:,9)  * dvr_mat(9,6)  &
                       +               &
               v(:,10) * dvr_mat(10,6)   

  v_scr(:,7) = v_scr(:,7)              &
                       +               &
               v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)  &
                       +               &
               v(:,8)  * dvr_mat(8,7)  &
                       +               &
               v(:,9)  * dvr_mat(9,7)  &
                       +               &
               v(:,10) * dvr_mat(10,7)   

  v_scr(:,8) = v_scr(:,8)              &
                       +               &
               v(:,1)  * dvr_mat(1,8)  &
                       +               &
               v(:,2)  * dvr_mat(2,8)  &
                       +               &
               v(:,3)  * dvr_mat(3,8)  &
                       +               &
               v(:,4)  * dvr_mat(4,8)  &
                       +               &
               v(:,5)  * dvr_mat(5,8)  &
                       +               &
               v(:,6)  * dvr_mat(6,8)  &
                       +               &
               v(:,7)  * dvr_mat(7,8)  &
                       +               &
               v(:,8)  * dvr_mat(8,8)  &
                       +               &
               v(:,9)  * dvr_mat(9,8)  &
                       +               &
               v(:,10) * dvr_mat(10,8)   

  v_scr(:,9) = v_scr(:,9)              &
                       +               &
               v(:,1)  * dvr_mat(1,9)  &
                       +               &
               v(:,2)  * dvr_mat(2,9)  &
                       +               &
               v(:,3)  * dvr_mat(3,9)  &
                       +               &
               v(:,4)  * dvr_mat(4,9)  &
                       +               &
               v(:,5)  * dvr_mat(5,9)  &
                       +               &
               v(:,6)  * dvr_mat(6,9)  &
                       +               &
               v(:,7)  * dvr_mat(7,9)  &
                       +               &
               v(:,8)  * dvr_mat(8,9)  &
                       +               &
               v(:,9)  * dvr_mat(9,9)  &
                       +               &
               v(:,10) * dvr_mat(10,9)   

  v_scr(:,10) = v_scr(:,10)            &
                       +               &
               v(:,1)  * dvr_mat(1,10)  &
                       +                &
               v(:,2)  * dvr_mat(2,10)  &
                       +                &
               v(:,3)  * dvr_mat(3,10)  &
                       +                &
               v(:,4)  * dvr_mat(4,10)  &
                       +                &
               v(:,5)  * dvr_mat(5,10)  &
                       +                &
               v(:,6)  * dvr_mat(6,10)  &
                       +                &
               v(:,7)  * dvr_mat(7,10)  &
                       +                &
               v(:,8)  * dvr_mat(8,10)  &
                       +                &
               v(:,9)  * dvr_mat(9,10)  &
                       +                &
               v(:,10) * dvr_mat(10,10)   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_10_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck v_v_mat_gen_z
!***begin prologue     v_v_mat_gen_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multipiles
!***                   for the general FEDVR HAmiltonian. 
!***                   The routine is not needed in 1D problems. 
!***                   In 2D, the i index runs over the number
!                      of points in the y coordinate while in 3D, it runs over
!                      the product of the number of points in the z and y coordinates.
!
!***references
!***routines called    ebcxx, ambcxx, apbcxx
!***end prologue   
!
  SUBROUTINE v_v_mat_gen_z(v,          &
                           v_scr,      &
                           dvr_mat,    & 
                           ni,nj,nk)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj, nk
  INTEGER                              :: j, k
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(nk,nk)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
!
!
  DO k=1,nk
     DO j=1,nk
        v_scr(:,j) = v_scr(:,j) + v(:,k) * dvr_mat(k,j) 
      END DO
   END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE v_v_mat_gen_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_2_z
!***begin prologue     v_v_mat_2_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagation matrix multiplies
!**                    for the special case of 2*2 matrices.
!
!***references
!***routines called
!***end prologue       v_v_mat_2_z
  SUBROUTINE v_v_mat_2_z(v,             &
                         v_scr,         &
                         dvr_mat,       & 
                         ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(2,2)               :: dvr_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  v_scr(:,1) = v_scr(:,1)                &
                       +                 &
               v(:,1)  * dvr_mat(1,1)    &
                       +                 &
               v(:,2)  * dvr_mat(2,1)    
 
  v_scr(:,2) = v_scr(:,2)                &
                       +                 &
               v(:,1)  * dvr_mat(1,2)    &
                       +                 &
               v(:,2)  * dvr_mat(2,2)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_2_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_3_z
!***begin prologue     v_v_mat_3_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_v_mat_3_z
!
  SUBROUTINE v_v_mat_3_z(v,          &
                         v_scr,      &
                         dvr_mat,    & 
                         ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(3,3)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       
  v_scr(:,1) = v_scr(:,1)                &
                       +                 &
               v(:,1)  * dvr_mat(1,1)    &
                       +                 &
               v(:,2)  * dvr_mat(2,1)    &
                       +                 &
               v(:,3)  * dvr_mat(3,1)    

  v_scr(:,2) = v_scr(:,2)                &
                       +                 &
               v(:,1)  * dvr_mat(1,2)    &
                       +                 &
               v(:,2)  * dvr_mat(2,2)    &
                       +                 &
               v(:,3)  * dvr_mat(3,2)    

  v_scr(:,3) = v_scr(:,3)                &
                       +                 &
               v(:,1)  * dvr_mat(1,3)    &
                       +                 &
               v(:,2)  * dvr_mat(2,3)    &
                       +                 &
               v(:,3)  * dvr_mat(3,3)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_3_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_4_z
!***begin prologue     v_v_mat_4_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_v_mat_4_z
!
  SUBROUTINE v_v_mat_4_z(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(4,4)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
  v_scr(:,1) = v_scr(:,1)                &
                       +                 &
               v(:,1)  * dvr_mat(1,1)    &
                       +                 &
               v(:,2)  * dvr_mat(2,1)    &
                       +                 &
               v(:,3)  * dvr_mat(3,1)    &
                       +                 &
               v(:,4)  * dvr_mat(4,1)    

  v_scr(:,2) = v_scr(:,2)                &
                       +                 &
               v(:,1)  * dvr_mat(1,2)    &
                       +                 &
               v(:,2)  * dvr_mat(2,2)    &
                       +                 &
               v(:,3)  * dvr_mat(3,2)    &
                       +                 &
               v(:,4)  * dvr_mat(4,2)    

  v_scr(:,3) = v_scr(:,3)                &
                       +                 &
               v(:,1)  * dvr_mat(1,3)    &
                       +                 &
               v(:,2)  * dvr_mat(2,3)    &
                       +                 &
               v(:,3)  * dvr_mat(3,3)    &
                       +                 &
               v(:,4)  * dvr_mat(4,3)    

  v_scr(:,4) = v_scr(:,4)                &
                       +                 &
               v(:,1)  * dvr_mat(1,4)    &
                       +                 &
               v(:,2)  * dvr_mat(2,4)    &
                       +                 &
               v(:,3)  * dvr_mat(3,4)    &
                       +                 &
               v(:,4)  * dvr_mat(4,4)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_4_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_5_z
!***begin prologue     v_v_mat_5_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_mat_5_z
!
  SUBROUTINE v_v_mat_5_z(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(5,5)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_5_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_6_z
!***begin prologue     v_v_mat_6_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_mat_6_z
!
  SUBROUTINE v_v_mat_6_z(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(6,6)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)    

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_6_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_7_z
!***begin prologue     v_v_mat_7_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_mat_7_z
!
  SUBROUTINE v_v_mat_7_z(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(7,7)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)    

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)    

  v_scr(:,7) = v_scr(:,7)              &
                       +               &
               v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_7_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_8_z
!***begin prologue     v_v_mat_8_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_mat_8_z
!
  SUBROUTINE v_v_mat_8_z(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(8,8)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)  &
                       +               &
               v(:,8)  * dvr_mat(8,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)  &
                       +               &
               v(:,8)  * dvr_mat(8,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)  &
                       +               &
               v(:,8)  * dvr_mat(8,5)    

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)  &
                       +               &
               v(:,8)  * dvr_mat(8,6)    

  v_scr(:,7) = v_scr(:,7)              &
                       +               &
               v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)  &
                       +               &
               v(:,8)  * dvr_mat(8,7)    

  v_scr(:,8) = v_scr(:,8)              &
                       +               &
               v(:,1)  * dvr_mat(1,8)  &
                       +               &
               v(:,2)  * dvr_mat(2,8)  &
                       +               &
               v(:,3)  * dvr_mat(3,8)  &
                       +               &
               v(:,4)  * dvr_mat(4,8)  &
                       +               &
               v(:,5)  * dvr_mat(5,8)  &
                       +               &
               v(:,6)  * dvr_mat(6,8)  &
                       +               &
               v(:,7)  * dvr_mat(7,8)  &
                       +               &
               v(:,8)  * dvr_mat(8,8)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_8_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_9_z
!***begin prologue     v_v_mat_9_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_mat_9_z
!
  SUBROUTINE v_v_mat_9_z(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(9,9)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)  &
                       +               &
               v(:,8)  * dvr_mat(8,1)  &
                       +               &
               v(:,9)  * dvr_mat(9,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)  &
                       +               &
               v(:,9)  * dvr_mat(9,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)  &
                       +               &
               v(:,8)  * dvr_mat(8,3)  &
                       +               &
               v(:,9)  * dvr_mat(9,3)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)  &
                       +               &
               v(:,8)  * dvr_mat(8,4)  &
                       +               &
               v(:,9)  * dvr_mat(9,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)  &
                       +               &
               v(:,8)  * dvr_mat(8,5)  &
                       +               &
               v(:,9)  * dvr_mat(9,5)    

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)  &
                       +               &
               v(:,8)  * dvr_mat(8,6)  &
                       +               &
               v(:,9)  * dvr_mat(9,6)    

  v_scr(:,7) = v_scr(:,7)              &
                       +               &
               v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)  &
                       +               &
               v(:,8)  * dvr_mat(8,7)  &
                       +               &
               v(:,9)  * dvr_mat(9,7)    

  v_scr(:,8) = v_scr(:,8)              &
                       +               &
               v(:,1)  * dvr_mat(1,8)  &
                       +               &
               v(:,2)  * dvr_mat(2,8)  &
                       +               &
               v(:,3)  * dvr_mat(3,8)  &
                       +               &
               v(:,4)  * dvr_mat(4,8)  &
                       +               &
               v(:,5)  * dvr_mat(5,8)  &
                       +               &
               v(:,6)  * dvr_mat(6,8)  &
                       +               &
               v(:,7)  * dvr_mat(7,8)  &
                       +               &
               v(:,8)  * dvr_mat(8,8)  &
                       +               &
               v(:,9)  * dvr_mat(9,8)    

  v_scr(:,9) = v_scr(:,9)              &
                       +               &
               v(:,1)  * dvr_mat(1,9)  &
                       +               &
               v(:,2)  * dvr_mat(2,9)  &
                       +               &
               v(:,3)  * dvr_mat(3,9)  &
                       +               &
               v(:,4)  * dvr_mat(4,9)  &
                       +               &
               v(:,5)  * dvr_mat(5,9)  &
                       +               &
               v(:,6)  * dvr_mat(6,9)  &
                       +               &
               v(:,7)  * dvr_mat(7,9)  &
                       +               &
               v(:,8)  * dvr_mat(8,9)  &
                       +               &
               v(:,9)  * dvr_mat(9,9)    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_9_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_10_z
!***begin prologue     v_v_mat_10_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_mat_10_z
!
  SUBROUTINE v_v_mat_10_z(v,        &
                          v_scr,    &
                          dvr_mat,  & 
                          ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  REAL*8, DIMENSION(10,10)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)  &
                       +               &
               v(:,8)  * dvr_mat(8,1)  &
                       +               &
               v(:,9)  * dvr_mat(9,1)  &
                       +               &
               v(:,10) * dvr_mat(10,1)   

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)  &
                       +               &
               v(:,9)  * dvr_mat(9,2)  &
                       +               &
               v(:,10) * dvr_mat(10,2)   

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)  &
                       +               &
               v(:,8)  * dvr_mat(8,3)  &
                       +               &
               v(:,9)  * dvr_mat(9,3)  &
                       +               &
               v(:,10) * dvr_mat(10,3)   

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)  &
                       +               &
               v(:,8)  * dvr_mat(8,4)  &
                       +               &
               v(:,9)  * dvr_mat(9,4)  &
                       +               &
               v(:,10) * dvr_mat(10,4)   

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)  &
                       +               &
               v(:,8)  * dvr_mat(8,5)  &
                       +               &
               v(:,9)  * dvr_mat(9,5)  &
                       +               &
               v(:,10) * dvr_mat(10,5)   

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)  &
                       +               &
               v(:,8)  * dvr_mat(8,6)  &
                       +               &
               v(:,9)  * dvr_mat(9,6)  &
                       +               &
               v(:,10) * dvr_mat(10,6)   

  v_scr(:,7) = v_scr(:,7)              &
                       +               &
               v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)  &
                       +               &
               v(:,8)  * dvr_mat(8,7)  &
                       +               &
               v(:,9)  * dvr_mat(9,7)  &
                       +               &
               v(:,10) * dvr_mat(10,7)   

  v_scr(:,8) = v_scr(:,8)              &
                       +               &
               v(:,1)  * dvr_mat(1,8)  &
                       +               &
               v(:,2)  * dvr_mat(2,8)  &
                       +               &
               v(:,3)  * dvr_mat(3,8)  &
                       +               &
               v(:,4)  * dvr_mat(4,8)  &
                       +               &
               v(:,5)  * dvr_mat(5,8)  &
                       +               &
               v(:,6)  * dvr_mat(6,8)  &
                       +               &
               v(:,7)  * dvr_mat(7,8)  &
                       +               &
               v(:,8)  * dvr_mat(8,8)  &
                       +               &
               v(:,9)  * dvr_mat(9,8)  &
                       +               &
               v(:,10) * dvr_mat(10,8)   

  v_scr(:,9) = v_scr(:,9)              &
                       +               &
               v(:,1)  * dvr_mat(1,9)  &
                       +               &
               v(:,2)  * dvr_mat(2,9)  &
                       +               &
               v(:,3)  * dvr_mat(3,9)  &
                       +               &
               v(:,4)  * dvr_mat(4,9)  &
                       +               &
               v(:,5)  * dvr_mat(5,9)  &
                       +               &
               v(:,6)  * dvr_mat(6,9)  &
                       +               &
               v(:,7)  * dvr_mat(7,9)  &
                       +               &
               v(:,8)  * dvr_mat(8,9)  &
                       +               &
               v(:,9)  * dvr_mat(9,9)  &
                       +               &
               v(:,10) * dvr_mat(10,9)   

  v_scr(:,10) = v_scr(:,10)            &
                       +               &
               v(:,1)  * dvr_mat(1,10)  &
                       +                &
               v(:,2)  * dvr_mat(2,10)  &
                       +                &
               v(:,3)  * dvr_mat(3,10)  &
                       +                &
               v(:,4)  * dvr_mat(4,10)  &
                       +                &
               v(:,5)  * dvr_mat(5,10)  &
                       +                &
               v(:,6)  * dvr_mat(6,10)  &
                       +                &
               v(:,7)  * dvr_mat(7,10)  &
                       +                &
               v(:,8)  * dvr_mat(8,10)  &
                       +                &
               v(:,9)  * dvr_mat(9,10)  &
                       +                &
               v(:,10) * dvr_mat(10,10)   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_10_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**begin prologue     v_out_so_mat_v_in
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       v_out_so_mat_v_in
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_gen_d
!***begin prologue     v_so_mat_v_gen_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propator matrix vector multiplies
!***                   for a general FEDVR Hamiltonian.
!***                   The parameter nj is a dummy and can take on values
!***                   consistent with any dimensional problem.  In a 1D
!***                   case nj=1, in 2D nj=nx and in 3D nj=ny*nx.
!
!***references
!***routines called    ebcxx, ambcxx, apbcxx
!***end prologue       
!
  SUBROUTINE v_so_mat_v_gen_d(v,             &
                              v_scr,         &
                              dvr_mat,       &
                              ni,nj,nk)
  IMPLICIT NONE
  INTEGER                              :: ni, nj, nk
  INTEGER                              :: i, k
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(nk,nk)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   v_scr(1:nk,:)=0.d0
   DO i=1,nk
      DO k=1,nk
         v_scr(k,:) = v_scr(k,:) + dvr_mat(k,i) * v(i,:) 
      END DO
   END DO
!
   v(1:nk,:) = v_scr(1:nk,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_gen_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!deck v_so_mat_v_2_d
!***begin prologue     v_so_mat_v_2_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies for
!***                   the special case of 2*2 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_2_d
!
  SUBROUTINE v_so_mat_v_2_d(v,         &
                            v_scr,     &
                            dvr_mat,   &
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(2,2)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   
  v_scr(2,:) = dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   
!
  v(1:2,:) = v_scr(1:2,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_2_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_3_d
!***begin prologue     v_so_mat_v_3_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_3_d
!
  SUBROUTINE v_so_mat_v_3_d(v,            &
                            v_scr,        &
                            dvr_mat,      &
                            ni,nj)
  IMPLICIT NONE
  INTEGER                                :: ni, nj
  REAL*8, DIMENSION(ni,nj)               :: v
  REAL*8, DIMENSION(ni,nj)               :: v_scr
  REAL*8, DIMENSION(3,3)                 :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   
  v_scr(2,:) = dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   
  v_scr(3,:) = dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   
!
  v(1:3,:) = v_scr(1:3,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_3_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_4_d
!***begin prologue     v_so_mat_v_4_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 4*4 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_4_d
!
  SUBROUTINE v_so_mat_v_4_d(v,             &
                            v_scr,         &
                            dvr_mat,       &
                            ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(4,4)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   &
                            +          &   
               dvr_mat(1,4) * v(4,:)   

  v_scr(2,:) = dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   &
                            +          &
               dvr_mat(2,4) * v(4,:)   

  v_scr(3,:) = dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   &
                            +          &
               dvr_mat(3,4) * v(4,:)   

  v_scr(4,:) = dvr_mat(4,1) * v(1,:)   &
                            +          &
               dvr_mat(4,2) * v(2,:)   &
                            +          &
               dvr_mat(4,3) * v(3,:)   &
                            +          &
               dvr_mat(4,4) * v(4,:)   
!
  v(1:4,:) = v_scr(1:4,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_4_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_5_d
!***begin prologue     v_so_mat_v_5_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_5_d
!
  SUBROUTINE v_so_mat_v_5_d(v,          &
                            v_scr,      &
                            dvr_mat,    &
                            ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(5,5)               :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   &
                            +          &   
               dvr_mat(1,4) * v(4,:)   &
                            +          &
               dvr_mat(1,5) * v(5,:)      
  v_scr(2,:) = dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   &
                            +          &
               dvr_mat(2,4) * v(4,:)   &
                            +          &
               dvr_mat(2,5) * v(5,:)   
  v_scr(3,:) = dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   &
                            +          &
               dvr_mat(3,4) * v(4,:)   &
                            +          &
               dvr_mat(3,5) * v(5,:)   
  v_scr(4,:) = dvr_mat(4,1) * v(1,:)   &
                            +          &
               dvr_mat(4,2) * v(2,:)   &
                            +          &
               dvr_mat(4,3) * v(3,:)   &
                            +          &
               dvr_mat(4,4) * v(4,:)   &
                            +          &
               dvr_mat(4,5) * v(5,:)    
  v_scr(5,:) = dvr_mat(5,1) * v(1,:)   &
                            +          &
               dvr_mat(5,2) * v(2,:)   &
                            +          &
               dvr_mat(5,3) * v(3,:)   &
                            +          &
               dvr_mat(5,4) * v(4,:)   &
                            +          &
               dvr_mat(5,5) * v(5,:)    
!
  v(1:5,:) = v_scr(1:5,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_5_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_6_d
!***begin prologue     v_so_mat_v_6_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_6_d
!
  SUBROUTINE v_so_mat_v_6_d(v,           &    
                            v_scr,       &
                            dvr_mat,     &
                            ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(6,6)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   &
                            +          &   
               dvr_mat(1,4) * v(4,:)   &
                            +          &
               dvr_mat(1,5) * v(5,:)   &   
                            +          &
               dvr_mat(1,6) * v(6,:)      
  v_scr(2,:) = dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   &
                            +          &   
               dvr_mat(2,4) * v(4,:)   &
                            +          &
               dvr_mat(2,5) * v(5,:)   &   
                            +          &
               dvr_mat(2,6) * v(6,:)      
  v_scr(3,:) = dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   &
                            +          &   
               dvr_mat(3,4) * v(4,:)   &
                            +          &
               dvr_mat(3,5) * v(5,:)   &   
                            +          &
               dvr_mat(3,6) * v(6,:)      
  v_scr(4,:) = dvr_mat(4,1) * v(1,:)   &
                            +          &
               dvr_mat(4,2) * v(2,:)   &
                            +          &
               dvr_mat(4,3) * v(3,:)   &
                            +          &   
               dvr_mat(4,4) * v(4,:)   &
                            +          &
               dvr_mat(4,5) * v(5,:)   &   
                            +          &
               dvr_mat(4,6) * v(6,:)      
  v_scr(5,:) = dvr_mat(5,1) * v(1,:)   &
                            +          &
               dvr_mat(5,2) * v(2,:)   &
                            +          &
               dvr_mat(5,3) * v(3,:)   &
                            +          &   
               dvr_mat(5,4) * v(4,:)   &
                            +          &
               dvr_mat(5,5) * v(5,:)   &   
                            +          &
               dvr_mat(5,6) * v(6,:)      

  v_scr(6,:) = dvr_mat(6,1) * v(1,:)   &
                            +          &
               dvr_mat(6,2) * v(2,:)   &
                            +          &
               dvr_mat(6,3) * v(3,:)   &
                            +          &   
               dvr_mat(6,4) * v(4,:)   &
                            +          &
               dvr_mat(6,5) * v(5,:)   &   
                            +          &
               dvr_mat(6,6) * v(6,:)      
!
  v(1:6,:) = v_scr(1:6,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_6_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_7_d
!***begin prologue     v_so_mat_v_7_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_7_d
!
  SUBROUTINE v_so_mat_v_7_d(v,           &    
                            v_scr,       &
                            dvr_mat,     &
                            ni,nj)
                     
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(7,7)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   &
                            +          &   
               dvr_mat(1,4) * v(4,:)   &
                            +          &
               dvr_mat(1,5) * v(5,:)   &   
                            +          &
               dvr_mat(1,6) * v(6,:)   &   
                            +          &
               dvr_mat(1,7) * v(7,:)      

  v_scr(2,:) = dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   &
                            +          &   
               dvr_mat(2,4) * v(4,:)   &
                            +          &
               dvr_mat(2,5) * v(5,:)   &   
                            +          &
               dvr_mat(2,6) * v(6,:)   &   
                            +          &
               dvr_mat(2,7) * v(7,:)      
  v_scr(3,:) = dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   &
                            +          &   
               dvr_mat(3,4) * v(4,:)   &
                            +          &
               dvr_mat(3,5) * v(5,:)   &   
                            +          &
               dvr_mat(3,6) * v(6,:)   &   
                            +          &
               dvr_mat(3,7) * v(7,:)      
  v_scr(4,:) = dvr_mat(4,1) * v(1,:)   &
                            +          &
               dvr_mat(4,2) * v(2,:)   &
                            +          &
               dvr_mat(4,3) * v(3,:)   &
                            +          &   
               dvr_mat(4,4) * v(4,:)   &
                            +          &
               dvr_mat(4,5) * v(5,:)   &   
                            +          &
               dvr_mat(4,6) * v(6,:)   &   
                            +          &
                dvr_mat(4,7) * v(7,:)      
  v_scr(5,:) = dvr_mat(5,1) * v(1,:)   &
                            +          &
               dvr_mat(5,2) * v(2,:)   &
                            +          &
               dvr_mat(5,3) * v(3,:)   &
                            +          &   
               dvr_mat(5,4) * v(4,:)   &
                            +          &
               dvr_mat(5,5) * v(5,:)   &   
                            +          &
               dvr_mat(5,6) * v(6,:)   &   
                            +          &
               dvr_mat(5,7) * v(7,:)      
  v_scr(6,:) = dvr_mat(6,1) * v(1,:)   &
                            +          &
               dvr_mat(6,2) * v(2,:)   &
                            +          &
               dvr_mat(6,3) * v(3,:)   &
                            +          &   
               dvr_mat(6,4) * v(4,:)   &
                            +          &
               dvr_mat(6,5) * v(5,:)   &   
                            +          &
               dvr_mat(6,6) * v(6,:)   &   
                            +          &
               dvr_mat(6,7) * v(7,:)      
  v_scr(7,:) = dvr_mat(7,1) * v(1,:)   &
                            +          &
               dvr_mat(7,2) * v(2,:)   &
                            +          &
               dvr_mat(7,3) * v(3,:)   &
                            +          &   
               dvr_mat(7,4) * v(4,:)   &
                            +          &
               dvr_mat(7,5) * v(5,:)   &   
                            +          &
               dvr_mat(7,6) * v(6,:)   &   
                            +          &
               dvr_mat(7,7) * v(7,:)      
!
  v(1:7,:) = v_scr(1:7,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_7_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_8_d
!***begin prologue     v_so_mat_v_8_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_8_d
!
  SUBROUTINE v_so_mat_v_8_d(v,             &
                            v_scr,         &
                            dvr_mat,       &
                            ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(8,8)               :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:)    = dvr_mat(1,1) * v(1,:)   &
                                 +        &
                  dvr_mat(1,2) * v(2,:)   &
                                 +        &
                  dvr_mat(1,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(1,4) * v(4,:)   &
                                 +        &
                  dvr_mat(1,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(1,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(1,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(1,8) * v(8,:)      

  v_scr(2,:) =    dvr_mat(2,1) * v(1,:)   &
                                 +        &
                  dvr_mat(2,2) * v(2,:)   &
                                 +        &
                  dvr_mat(2,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(2,4) * v(4,:)   &
                                 +        &
                  dvr_mat(2,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(2,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(2,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(2,8) * v(8,:)      

  v_scr(3,:) =    dvr_mat(3,1) * v(1,:)   &
                                 +        &
                  dvr_mat(3,2) * v(2,:)   &
                                 +        &
                  dvr_mat(3,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(3,4) * v(4,:)   &
                                 +        &
                  dvr_mat(3,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(3,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(3,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(3,8) * v(8,:)      

  v_scr(4,:) =    dvr_mat(4,1) * v(1,:)   &
                                 +        &
                  dvr_mat(4,2) * v(2,:)   &
                                 +        &
                  dvr_mat(4,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(4,4) * v(4,:)   &
                                 +        &
                  dvr_mat(4,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(4,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(4,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(4,8) * v(8,:)      

  v_scr(5,:) =    dvr_mat(5,1) * v(1,:)   &
                                 +        &
                   dvr_mat(5,2) * v(2,:)  &
                                 +        &
                  dvr_mat(5,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(5,4) * v(4,:)   &
                                 +        &
                  dvr_mat(5,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(5,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(5,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(5,8) * v(8,:)     
   
  v_scr(6,:) =    dvr_mat(6,1) * v(1,:)   &
                                 +        &
                  dvr_mat(6,2) * v(2,:)   &
                                 +        &
                  dvr_mat(6,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(6,4) * v(4,:)   &
                                 +        &
                  dvr_mat(6,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(6,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(6,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(6,8) * v(8,:)        

  v_scr(7,:) =    dvr_mat(7,1) * v(1,:)   &
                                 +        &
                  dvr_mat(7,2) * v(2,:)   &
                                 +        &
                  dvr_mat(7,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(7,4) * v(4,:)   &
                                 +        &
                  dvr_mat(7,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(7,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(7,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(7,8) * v(8,:)        

  v_scr(8,:) =    dvr_mat(8,1) * v(1,:)   &
                                 +        &
                  dvr_mat(8,2) * v(2,:)   &
                                 +        &
                  dvr_mat(8,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(8,4) * v(4,:)   &
                                 +        &
                  dvr_mat(8,5) * v(5,:)   &  
                                 +        &
                  dvr_mat(8,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(8,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(8,8) * v(8,:)        

!

  v(1:8,:) = v_scr(1:8,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_8_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_9_d
!***begin prologue     v_so_mat_v_9_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_9_d
!
  SUBROUTINE v_so_mat_v_9_d(v,           &
                            v_scr,       &
                            dvr_mat,     &
                            ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(9,9)               :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) =    dvr_mat(1,1) * v(1,:)   &
                                 +        &
                  dvr_mat(1,2) * v(2,:)   &
                                 +        &
                  dvr_mat(1,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(1,4) * v(4,:)   &
                                 +        &
                  dvr_mat(1,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(1,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(1,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(1,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(1,9) * v(9,:)      
  v_scr(2,:) =    dvr_mat(2,1) * v(1,:)   &
                                 +        &
                  dvr_mat(2,2) * v(2,:)   &
                                 +        &
                  dvr_mat(2,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(2,4) * v(4,:)   &
                                 +        &
                  dvr_mat(2,5) * v(5,:)   &
                                 +        &
                  dvr_mat(2,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(2,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(2,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(2,9) * v(9,:)      

  v_scr(3,:) =    dvr_mat(3,1) * v(1,:)   &
                                 +        &
                  dvr_mat(3,2) * v(2,:)   &
                                 +        &
                  dvr_mat(3,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(3,4) * v(4,:)   &
                                 +        &
                  dvr_mat(3,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(3,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(3,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(3,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(3,9) * v(9,:)      

  v_scr(4,:) =    dvr_mat(4,1) * v(1,:)   &
                                 +        &
                  dvr_mat(4,2) * v(2,:)   &
                                 +        &
                  dvr_mat(4,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(4,4) * v(4,:)   &
                                 +        &
                  dvr_mat(4,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(4,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(4,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(4,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(4,9) * v(9,:)      

  v_scr(5,:) =    dvr_mat(5,1) * v(1,:)   &
                                 +        &
                  dvr_mat(5,2) * v(2,:)   &
                                 +        &
                  dvr_mat(5,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(5,4) * v(4,:)   &
                                 +        &
                  dvr_mat(5,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(5,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(5,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(5,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(5,9) * v(9,:)      

  v_scr(6,:) =    dvr_mat(6,1) * v(1,:)   &
                                 +        &
                  dvr_mat(6,2) * v(2,:)   &
                                 +        &
                  dvr_mat(6,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(6,4) * v(4,:)   &
                                 +        &
                  dvr_mat(6,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(6,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(6,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(6,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(6,9) * v(9,:)      

  v_scr(7,:) =    dvr_mat(7,1) * v(1,:)   &
                                 +        &
                  dvr_mat(7,2) * v(2,:)   &
                                 +        &
                  dvr_mat(7,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(7,4) * v(4,:)   &
                                 +        &
                  dvr_mat(7,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(7,6) * v(6,:)   &  
                                 +        &
                  dvr_mat(7,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(7,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(7,9) * v(9,:)      

  v_scr(8,:) =    dvr_mat(8,1) * v(1,:)   &
                                 +        &
                  dvr_mat(8,2) * v(2,:)   &
                                 +        &
                  dvr_mat(8,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(8,4) * v(4,:)   &
                                 +        &
                  dvr_mat(8,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(8,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(8,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(8,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(8,9) * v(9,:)      

  v_scr(9,:) =    dvr_mat(9,1) * v(1,:)   &
                                 +        &
                  dvr_mat(9,2) * v(2,:)   &
                                 +        &
                  dvr_mat(9,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(9,4) * v(4,:)   &
                                 +        &
                  dvr_mat(9,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(9,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(9,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(9,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(9,9) * v(9,:)      
!
  v(1:9,:) = v_scr(1:9,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_9_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_10_d
!***begin prologue     v_so_mat_v_10_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_10_d
!
  SUBROUTINE v_so_mat_v_10_d(v,          &
                             v_scr,      &
                             dvr_mat,    &
                             ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(10,10)             :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) =    dvr_mat(1,1) * v(1,:)   &
                                 +        &
                  dvr_mat(1,2) * v(2,:)   &
                                 +        &
                  dvr_mat(1,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(1,4) * v(4,:)   &
                                 +        &
                  dvr_mat(1,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(1,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(1,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(1,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(1,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(1,10) * v(10,:)    

  v_scr(2,:) =    dvr_mat(2,1) * v(1,:)   &
                                 +        &
                  dvr_mat(2,2) * v(2,:)   &
                                 +        &
                  dvr_mat(2,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(2,4) * v(4,:)   &
                                 +        &
                  dvr_mat(2,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(2,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(2,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(2,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(2,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(2,10) * v(10,:)    

  v_scr(3,:) =    dvr_mat(3,1) * v(1,:)   &
                                 +        &
                  dvr_mat(3,2) * v(2,:)   &
                                 +        &
                  dvr_mat(3,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(3,4) * v(4,:)   &
                                 +        &
                  dvr_mat(3,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(3,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(3,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(3,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(3,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(3,10) * v(10,:)    

  v_scr(4,:) =    dvr_mat(4,1) * v(1,:)   &
                                 +        &
                  dvr_mat(4,2) * v(2,:)   &
                                 +        &
                  dvr_mat(4,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(4,4) * v(4,:)   &
                                 +        &
                  dvr_mat(4,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(4,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(4,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(4,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(4,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(4,10) * v(10,:)    

  v_scr(5,:) =    dvr_mat(5,1) * v(1,:)   &
                                 +        &
                  dvr_mat(5,2) * v(2,:)  &
                                 +        &
                  dvr_mat(5,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(5,4) * v(4,:)   &
                                 +        &
                  dvr_mat(5,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(5,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(5,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(5,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(5,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(5,10) * v(10,:)    
  v_scr(6,:) =    dvr_mat(6,1) * v(1,:)   &
                                 +        &
                  dvr_mat(6,2) * v(2,:)   &
                                 +        &
                  dvr_mat(6,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(6,4) * v(4,:)   &
                                 +        &
                  dvr_mat(6,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(6,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(6,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(6,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(6,9) * v(9,:)   &  
                                 +        &
                  dvr_mat(6,10) * v(10,:)     

  v_scr(7,:) =    dvr_mat(7,1) * v(1,:)   &
                                 +        &
                  dvr_mat(7,2) * v(2,:)   &
                                 +        &
                  dvr_mat(7,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(7,4) * v(4,:)   &
                                 +        &
                  dvr_mat(7,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(7,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(7,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(7,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(7,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(7,10) * v(10,:)   

  v_scr(8,:) =    dvr_mat(8,1) * v(1,:)   &
                                 +        &
                  dvr_mat(8,2) * v(2,:)   &
                                 +        &
                  dvr_mat(8,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(8,4) * v(4,:)   &
                                 +        &
                  dvr_mat(8,5) * v(5,:)   &  
                                 +        &
                  dvr_mat(8,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(8,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(8,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(8,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(8,10) * v(10,:)      

  v_scr(9,:) =    dvr_mat(9,1) * v(1,:)   &
                                 +        &
                  dvr_mat(9,2) * v(2,:)   &
                                 +        &
                  dvr_mat(9,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(9,4) * v(4,:)   &
                                 +        &
                  dvr_mat(9,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(9,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(9,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(9,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(9,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(9,10) * v(10,:)    

  v_scr(10,:) =   dvr_mat(10,1) * v(1,:) &
                                   +     &
                  dvr_mat(10,2) * v(2,:) &
                                   +     &
                  dvr_mat(10,3) * v(3,:) &
                                   +     &   
                  dvr_mat(10,4) * v(4,:) &
                                   +     &
                  dvr_mat(10,5) * v(5,:) &   
                                   +     &
                  dvr_mat(10,6) * v(6,:) &   
                                   +     &
                  dvr_mat(10,7) * v(7,:) &   
                                   +     &
                  dvr_mat(10,8) * v(8,:) &   
                                   +     &
                  dvr_mat(10,9) * v(9,:) &   
                                   +     &
                  dvr_mat(10,10) * v(10,:)    
!
  v(1:10,:) = v_scr(1:10,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_10_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_11_d
!***begin prologue     v_so_mat_v_11_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_11_d
!
  SUBROUTINE v_so_mat_v_11_d(v,          &
                             v_scr,      &
                             dvr_mat,    &
                             ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(11,11)             :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) =    dvr_mat(1,1) * v(1,:)   &
                                 +        &
                  dvr_mat(1,2) * v(2,:)   &
                                 +        &
                  dvr_mat(1,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(1,4) * v(4,:)   &
                                 +        &
                  dvr_mat(1,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(1,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(1,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(1,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(1,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(1,10) * v(10,:) &   
                                 +        &
                  dvr_mat(1,11) * v(11,:)    

  v_scr(2,:) =    dvr_mat(2,1) * v(1,:)   &
                                 +        &
                  dvr_mat(2,2) * v(2,:)   &
                                 +        &
                  dvr_mat(2,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(2,4) * v(4,:)   &
                                 +        &
                  dvr_mat(2,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(2,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(2,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(2,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(2,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(2,10) * v(10,:) &   
                                 +        &
                  dvr_mat(2,11) * v(11,:)    

  v_scr(3,:) =    dvr_mat(3,1) * v(1,:)   &
                                 +        &
                  dvr_mat(3,2) * v(2,:)   &
                                 +        &
                  dvr_mat(3,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(3,4) * v(4,:)   &
                                 +        &
                  dvr_mat(3,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(3,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(3,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(3,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(3,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(3,10) * v(10,:) &   
                                 +        &
                  dvr_mat(3,11) * v(11,:)    

  v_scr(4,:) =    dvr_mat(4,1) * v(1,:)   &
                                 +        &
                  dvr_mat(4,2) * v(2,:)   &
                                 +        &
                  dvr_mat(4,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(4,4) * v(4,:)   &
                                 +        &
                  dvr_mat(4,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(4,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(4,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(4,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(4,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(4,10) * v(10,:) &   
                                 +        &
                  dvr_mat(4,11) * v(11,:)    

  v_scr(5,:) =    dvr_mat(5,1) * v(1,:)   &
                                 +        &
                  dvr_mat(5,2) * v(2,:)  &
                                 +        &
                  dvr_mat(5,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(5,4) * v(4,:)   &
                                 +        &
                  dvr_mat(5,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(5,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(5,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(5,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(5,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(5,10) * v(10,:) &   
                                 +        &
                  dvr_mat(5,11) * v(11,:)    

  v_scr(6,:) =    dvr_mat(6,1) * v(1,:)   &
                                 +        &
                  dvr_mat(6,2) * v(2,:)   &
                                 +        &
                  dvr_mat(6,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(6,4) * v(4,:)   &
                                 +        &
                  dvr_mat(6,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(6,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(6,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(6,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(6,9) * v(9,:)   &  
                                 +        &
                  dvr_mat(6,10) * v(10,:) &     
                                 +        &
                  dvr_mat(6,11) * v(11,:)    

  v_scr(7,:) =    dvr_mat(7,1) * v(1,:)   &
                                 +        &
                  dvr_mat(7,2) * v(2,:)   &
                                 +        &
                  dvr_mat(7,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(7,4) * v(4,:)   &
                                 +        &
                  dvr_mat(7,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(7,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(7,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(7,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(7,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(7,10) * v(10,:) &  
                                 +        &
                  dvr_mat(7,11) * v(11,:)    

  v_scr(8,:) =    dvr_mat(8,1) * v(1,:)   &
                                 +        &
                  dvr_mat(8,2) * v(2,:)   &
                                 +        &
                  dvr_mat(8,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(8,4) * v(4,:)   &
                                 +        &
                  dvr_mat(8,5) * v(5,:)   &  
                                 +        &
                  dvr_mat(8,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(8,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(8,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(8,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(8,10) * v(10,:) &     
                                 +        &
                  dvr_mat(9,11) * v(11,:)    

  v_scr(9,:) =    dvr_mat(9,1) * v(1,:)   &
                                 +        &
                  dvr_mat(9,2) * v(2,:)   &
                                 +        &
                  dvr_mat(9,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(9,4) * v(4,:)   &
                                 +        &
                  dvr_mat(9,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(9,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(9,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(9,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(9,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(9,10) * v(10,:) &   
                                 +        &
                  dvr_mat(9,11) * v(11,:)    

  v_scr(10,:) =   dvr_mat(10,1) * v(1,:) &
                                   +     &
                  dvr_mat(10,2) * v(2,:) &
                                   +     &
                  dvr_mat(10,3) * v(3,:) &
                                   +     &   
                  dvr_mat(10,4) * v(4,:) &
                                   +     &
                  dvr_mat(10,5) * v(5,:) &   
                                   +     &
                  dvr_mat(10,6) * v(6,:) &   
                                   +     &
                  dvr_mat(10,7) * v(7,:) &   
                                   +     &
                  dvr_mat(10,8) * v(8,:) &   
                                   +     &
                  dvr_mat(10,9) * v(9,:) &   
                                   +     &
                  dvr_mat(10,10) * v(10,:)&   
                                 +        &
                  dvr_mat(10,11) * v(11,:)     

  v_scr(11,:) =   dvr_mat(11,1) * v(1,:) &
                                   +     &
                  dvr_mat(11,2) * v(2,:) &
                                   +     &
                  dvr_mat(11,3) * v(3,:) &
                                   +     &   
                  dvr_mat(11,4) * v(4,:) &
                                   +     &
                  dvr_mat(11,5) * v(5,:) &   
                                   +     &
                  dvr_mat(11,6) * v(6,:) &   
                                   +     &
                  dvr_mat(11,7) * v(7,:) &   
                                   +     &
                  dvr_mat(11,8) * v(8,:) &   
                                   +     &
                  dvr_mat(11,9) * v(9,:) &   
                                   +     &
                  dvr_mat(11,10) * v(10,:)&   
                                 +        &
                  dvr_mat(11,11) * v(11,:)     
!
  v(1:11,:) = v_scr(1:11,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_11_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_gen_z
!***begin prologue     v_so_mat_v_gen_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propator matrix vector multiplies
!***                   for a general FEDVR Hamiltonian.
!***                   The parameter nj is a dummy and can take on values
!***                   consistent with any dimensional problem.  In a 1D
!***                   case nj=1, in 2D nj=nx and in 3D nj=ny*nx.
!
!***references
!***routines called    ebcxx, ambcxx, apbcxx
!***end prologue       
!
  SUBROUTINE v_so_mat_v_gen_z(v,                                   &
                              v_scr,                               &
                              dvr_mat,                             &
                              ni,nj,nk)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, nk
  INTEGER                                  :: i, k
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(nk,nk)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  v_scr(1:nk,:)= (0.d0,0.d0)
  DO i=1,nk
     DO k=1,nk
        v_scr(k,:) = v_scr(k,:) + dvr_mat(k,i) * v(i,:)     
     END DO
  END DO
!
  v(1:nk,:) = v_scr(1:nk,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE  v_so_mat_v_gen_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_2_z
!***begin prologue     v_so_mat_v_2_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies for
!***                   the special case of 2*2 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_2_z
!
  SUBROUTINE v_so_mat_v_2_z(v,          &
                            v_scr,      &
                            dvr_mat,    &
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(2,2)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)     &
                            +            &
               dvr_mat(1,2) * v(2,:)  
  v_scr(2,:) = dvr_mat(2,1) * v(1,:)     &
                            +            &
               dvr_mat(2,2) * v(2,:)     
!
  v(1:2,:) = v_scr(1:2,:)
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_2_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_3_z
!***begin prologue     v_so_mat_v_3_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_3_z
!
  SUBROUTINE v_so_mat_v_3_z(v,          &
                            v_scr,      &
                            dvr_mat,    &
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(3,3)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)     &
                            +            &
               dvr_mat(1,2) * v(2,:)     &
                            +            &
               dvr_mat(1,3) * v(3,:)     
  v_scr(2,:) = dvr_mat(2,1) * v(1,:)     &
                            +            &
               dvr_mat(2,2) * v(2,:)     &
                            +            &
               dvr_mat(2,3) * v(3,:)     
  v_scr(3,:) = dvr_mat(3,1) * v(1,:)     &
                            +            &
               dvr_mat(3,2) * v(2,:)     &
                            +            &
               dvr_mat(3,3) * v(3,:)     
!
  v(1:3,:) = v_scr(1:3,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_3_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_4_z
!***begin prologue     v_so_mat_v_4_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_4_z
!
  SUBROUTINE v_so_mat_v_4_z(v,          &
                            v_scr,      &
                            dvr_mat,    &
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(4,4)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)     &
                            +            &
               dvr_mat(1,2) * v(2,:)     &
                            +            &
               dvr_mat(1,3) * v(3,:)     &
                            +            &
               dvr_mat(1,4) * v(4,:)     

  v_scr(2,:) = dvr_mat(2,1) * v(1,:)     &
                            +            &
               dvr_mat(2,2) * v(2,:)     &
                            +            &
               dvr_mat(2,3) * v(3,:)     &
                            +            &
               dvr_mat(2,4) * v(4,:)     

  v_scr(3,:) = dvr_mat(3,1) * v(1,:)     &
                            +            &
               dvr_mat(3,2) * v(2,:)     &
                            +            &
               dvr_mat(3,3) * v(3,:)     &
                            +            &
               dvr_mat(3,4) * v(4,:)     

  v_scr(4,:) = dvr_mat(4,1) * v(1,:)     &
                            +            &
               dvr_mat(4,2) * v(2,:)     &
                            +            &
               dvr_mat(4,3) * v(3,:)     &
                            +            &
               dvr_mat(4,4) * v(4,:)     
!
  v(1:4,:) = v_scr(1:4,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_4_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_5_z
!***begin prologue     v_so_mat_v_5_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_5_z
!
  SUBROUTINE v_so_mat_v_5_z(v,          &
                            v_scr,      &
                            dvr_mat,    &
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(5,5)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)     &
                            +            &
               dvr_mat(1,2) * v(2,:)     &
                            +            &
               dvr_mat(1,3) * v(3,:)     &
                            +            &
               dvr_mat(1,4) * v(4,:)     &
                            +            &
               dvr_mat(1,5) * v(5,:)     

  v_scr(2,:) = dvr_mat(2,1) * v(1,:)     &
                            +            &
               dvr_mat(2,2) * v(2,:)     &
                            +            &
               dvr_mat(2,3) * v(3,:)     &
                            +            &
               dvr_mat(2,4) * v(4,:)     &
                            +            &
               dvr_mat(2,5) * v(5,:)     

  v_scr(3,:) = dvr_mat(3,1) * v(1,:)     &
                            +            &
               dvr_mat(3,2) * v(2,:)     &
                            +            &
               dvr_mat(3,3) * v(3,:)     &
                            +            &
               dvr_mat(3,4) * v(4,:)     &
                            +            &
               dvr_mat(3,5) * v(5,:)     

  v_scr(4,:) = dvr_mat(4,1) * v(1,:)     &
                            +            &
               dvr_mat(4,2) * v(2,:)     &
                            +            &
               dvr_mat(4,3) * v(3,:)     &
                            +            &
               dvr_mat(4,4) * v(4,:)     &
                            +            &
               dvr_mat(4,5) * v(5,:)     

  v_scr(5,:) = dvr_mat(5,1) * v(1,:)     &
                            +            &
               dvr_mat(5,2) * v(2,:)     &
                            +            &
               dvr_mat(5,3) * v(3,:)     &
                            +            &
               dvr_mat(5,4) * v(4,:)     &
                            +            &
               dvr_mat(5,5) * v(5,:)     

!
  v(1:5,:) = v_scr(1:5,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_5_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_6_z
!***begin prologue     v_so_mat_v_6_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_6_z
!
  SUBROUTINE v_so_mat_v_6_z(v,          &
                            v_scr,      &
                            dvr_mat,    &
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  COMPLEX*16, DIMENSION(6,6)           :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)     &
                            +            &
               dvr_mat(1,2) * v(2,:)     &
                            +            &
               dvr_mat(1,3) * v(3,:)     &
                            +            &
               dvr_mat(1,4) * v(4,:)     &
                            +            &
               dvr_mat(1,5) * v(5,:)     &
                            +            &
               dvr_mat(1,6) * v(6,:)     

  v_scr(2,:) = dvr_mat(2,1) * v(1,:)     &
                            +            &
               dvr_mat(2,2) * v(2,:)     &
                            +            &
               dvr_mat(2,3) * v(3,:)     &
                            +            &
               dvr_mat(2,4) * v(4,:)     &
                            +            &
               dvr_mat(2,5) * v(5,:)     &
                            +            &
               dvr_mat(2,6) * v(6,:)     

  v_scr(3,:) = dvr_mat(3,1) * v(1,:)     &
                            +            &
               dvr_mat(3,2) * v(2,:)     &
                            +            &
               dvr_mat(3,3) * v(3,:)     &
                            +            &
               dvr_mat(3,4) * v(4,:)     &
                            +            &
               dvr_mat(3,5) * v(5,:)     &
                            +            &
               dvr_mat(3,6) * v(6,:)     

  v_scr(4,:) = dvr_mat(4,1) * v(1,:)     &
                            +            &
               dvr_mat(4,2) * v(2,:)     &
                            +            &
               dvr_mat(4,3) * v(3,:)     &
                            +            &
               dvr_mat(4,4) * v(4,:)     &
                            +            &
               dvr_mat(4,5) * v(5,:)     &
                            +            &
               dvr_mat(4,6) * v(6,:)     

  v_scr(5,:) = dvr_mat(5,1) * v(1,:)     &
                            +            &
               dvr_mat(5,2) * v(2,:)     &
                            +            &
               dvr_mat(5,3) * v(3,:)     &
                            +            &
               dvr_mat(5,4) * v(4,:)     &
                            +            &
               dvr_mat(5,5) * v(5,:)     &
                            +            &
               dvr_mat(5,6) * v(6,:)     

  v_scr(6,:) = dvr_mat(6,1) * v(1,:)     &
                            +            &
               dvr_mat(6,2) * v(2,:)     &
                            +            &
               dvr_mat(6,3) * v(3,:)     &
                            +            &
               dvr_mat(6,4) * v(4,:)     &
                            +            &
               dvr_mat(6,5) * v(5,:)     &
                            +            &
               dvr_mat(6,6) * v(6,:)     

!
  v(1:6,:) = v_scr(1:6,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_6_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_7_z
!***begin prologue     v_so_mat_v_7_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_7_z
!
  SUBROUTINE v_so_mat_v_7_z(v,          &
                            v_scr,      &
                            dvr_mat,    &
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(7,7)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)     &
                            +            &
               dvr_mat(1,2) * v(2,:)     &
                            +            &
               dvr_mat(1,3) * v(3,:)     &
                            +            &
               dvr_mat(1,4) * v(4,:)     &
                            +            &
               dvr_mat(1,5) * v(5,:)     &
                            +            &
               dvr_mat(1,6) * v(6,:)     &
                            +            &
               dvr_mat(1,7) * v(7,:)     

  v_scr(2,:) = dvr_mat(2,1) * v(1,:)     &
                            +            &
               dvr_mat(2,2) * v(2,:)     &
                            +            &
               dvr_mat(2,3) * v(3,:)     &
                            +            &
               dvr_mat(2,4) * v(4,:)     &
                            +            &
               dvr_mat(2,5) * v(5,:)     &
                            +            &
               dvr_mat(2,6) * v(6,:)     &
                            +            &
               dvr_mat(2,7) * v(7,:)     

  v_scr(3,:) = dvr_mat(3,1) * v(1,:)     &
                            +            &
               dvr_mat(3,2) * v(2,:)     &
                            +            &
               dvr_mat(3,3) * v(3,:)     &
                            +            &
               dvr_mat(3,4) * v(4,:)     &
                            +            &
               dvr_mat(3,5) * v(5,:)     &
                            +            &
               dvr_mat(3,6) * v(6,:)     &
                            +            &
               dvr_mat(3,7) * v(7,:)     

  v_scr(4,:) = dvr_mat(4,1) * v(1,:)     &
                            +            &
               dvr_mat(4,2) * v(2,:)     &
                            +            &
               dvr_mat(4,3) * v(3,:)     &
                            +            &
               dvr_mat(4,4) * v(4,:)     &
                            +            &
               dvr_mat(4,5) * v(5,:)     &
                            +            &
               dvr_mat(4,6) * v(6,:)     &
                            +            &
               dvr_mat(4,7) * v(7,:)     

  v_scr(5,:) = dvr_mat(5,1) * v(1,:)     &
                            +            &
               dvr_mat(5,2) * v(2,:)     &
                            +            &
               dvr_mat(5,3) * v(3,:)     &
                            +            &
               dvr_mat(5,4) * v(4,:)     &
                            +            &
               dvr_mat(5,5) * v(5,:)     &
                            +            &
               dvr_mat(5,6) * v(6,:)     &
                            +            &
               dvr_mat(5,7) * v(7,:)     

  v_scr(6,:) = dvr_mat(6,1) * v(1,:)     &
                            +            &
               dvr_mat(6,2) * v(2,:)     &
                            +            &
               dvr_mat(6,3) * v(3,:)     &
                            +            &
               dvr_mat(6,4) * v(4,:)     &
                            +            &
               dvr_mat(6,5) * v(5,:)     &
                            +            &
               dvr_mat(6,6) * v(6,:)     &
                            +            &
               dvr_mat(6,7) * v(7,:)     

  v_scr(7,:) = dvr_mat(7,1) * v(1,:)     &
                            +            &
               dvr_mat(7,2) * v(2,:)     &
                            +            &
               dvr_mat(7,3) * v(3,:)     &
                            +            &
               dvr_mat(7,4) * v(4,:)     &
                            +            &
               dvr_mat(7,5) * v(5,:)     &
                            +            &
               dvr_mat(7,6) * v(6,:)     &
                            +            &
               dvr_mat(7,7) * v(7,:)     

!
  v(1:7,:) = v_scr(1:7,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_7_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_8_z
!***begin prologue     v_so_mat_v_8_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_8_z
!
  SUBROUTINE v_so_mat_v_8_z(v,               &
                            v_scr,           &
                            dvr_mat,         &
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(8,8)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)     &
                            +            &
               dvr_mat(1,2) * v(2,:)     &
                            +            &
               dvr_mat(1,3) * v(3,:)     &
                            +            &
               dvr_mat(1,4) * v(4,:)     &
                            +            &
               dvr_mat(1,5) * v(5,:)     &
                            +            &
               dvr_mat(1,6) * v(6,:)     &
                            +            &
               dvr_mat(1,7) * v(7,:)     &
                            +            &
               dvr_mat(1,8) * v(8,:)     

  v_scr(2,:) = dvr_mat(2,1) * v(1,:)     &
                            +            &
               dvr_mat(2,2) * v(2,:)     &
                            +            &
               dvr_mat(2,3) * v(3,:)     &
                            +            &
               dvr_mat(2,4) * v(4,:)     &
                            +            &
               dvr_mat(2,5) * v(5,:)     &
                            +            &
               dvr_mat(2,6) * v(6,:)     &
                            +            &
               dvr_mat(2,7) * v(7,:)     &
                            +            &
               dvr_mat(2,8) * v(8,:)     

  v_scr(3,:) = dvr_mat(3,1) * v(1,:)     &
                            +            &
               dvr_mat(3,2) * v(2,:)     &
                            +            &
               dvr_mat(3,3) * v(3,:)     &
                            +            &
               dvr_mat(3,4) * v(4,:)     &
                            +            &
               dvr_mat(3,5) * v(5,:)     &
                            +            &
               dvr_mat(3,6) * v(6,:)     &
                            +            &
               dvr_mat(3,7) * v(7,:)     &
                            +            &
               dvr_mat(3,8) * v(8,:)     

  v_scr(4,:) = dvr_mat(4,1) * v(1,:)     &
                            +            &
               dvr_mat(4,2) * v(2,:)     &
                            +            &
               dvr_mat(4,3) * v(3,:)     &
                            +            &
               dvr_mat(4,4) * v(4,:)     &
                            +            &
               dvr_mat(4,5) * v(5,:)     &
                            +            &
               dvr_mat(4,6) * v(6,:)     &
                            +            &
               dvr_mat(4,7) * v(7,:)     &
                            +            &
               dvr_mat(4,8) * v(8,:)     

  v_scr(5,:) = dvr_mat(5,1) * v(1,:)     &
                            +            &
               dvr_mat(5,2) * v(2,:)     &
                            +            &
               dvr_mat(5,3) * v(3,:)     &
                            +            &
               dvr_mat(5,4) * v(4,:)     &
                            +            &
               dvr_mat(5,5) * v(5,:)     &
                            +            &
               dvr_mat(5,6) * v(6,:)     &
                            +            &
               dvr_mat(5,7) * v(7,:)     &
                            +            &
               dvr_mat(5,8) * v(8,:)     

  v_scr(6,:) = dvr_mat(6,1) * v(1,:)     &
                            +            &
               dvr_mat(6,2) * v(2,:)     &
                            +            &
               dvr_mat(6,3) * v(3,:)     &
                            +            &
               dvr_mat(6,4) * v(4,:)     &
                            +            &
               dvr_mat(6,5) * v(5,:)     &
                            +            &
               dvr_mat(6,6) * v(6,:)     &
                            +            &
               dvr_mat(6,7) * v(7,:)     &
                            +            &
               dvr_mat(6,8) * v(8,:)     

  v_scr(7,:) = dvr_mat(7,1) * v(1,:)     &
                            +            &
               dvr_mat(7,2) * v(2,:)     &
                            +            &
               dvr_mat(7,3) * v(3,:)     &
                            +            &
               dvr_mat(7,4) * v(4,:)     &
                            +            &
               dvr_mat(7,5) * v(5,:)     &
                            +            &
               dvr_mat(7,6) * v(6,:)     &
                            +            &
               dvr_mat(7,7) * v(7,:)     &
                            +            &
               dvr_mat(7,8) * v(8,:)     

  v_scr(8,:) = dvr_mat(8,1) * v(1,:)     &
                            +            &
               dvr_mat(8,2) * v(2,:)     &
                            +            &
               dvr_mat(8,3) * v(3,:)     &
                            +            &
               dvr_mat(8,4) * v(4,:)     &
                            +            &
               dvr_mat(8,5) * v(5,:)     &
                            +            &
               dvr_mat(8,6) * v(6,:)     &
                            +            &
               dvr_mat(8,7) * v(7,:)     &
                            +            &
               dvr_mat(8,8) * v(8,:)     

!
  v(1:8,:) = v_scr(1:8,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_8_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_9_z
!***begin prologue     v_so_mat_v_9_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_9_z
!
  SUBROUTINE v_so_mat_v_9_z(v,               &
                            v_scr,           &
                            dvr_mat,         &
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(9,9)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)     &
                            +            &
               dvr_mat(1,2) * v(2,:)     &
                            +            &
               dvr_mat(1,3) * v(3,:)     &
                            +            &
               dvr_mat(1,4) * v(4,:)     &
                            +            &
               dvr_mat(1,5) * v(5,:)     &
                            +            &
               dvr_mat(1,6) * v(6,:)     &
                            +            &
               dvr_mat(1,7) * v(7,:)     &
                            +            &
               dvr_mat(1,8) * v(8,:)     &
                            +            &
               dvr_mat(1,9) * v(9,:)     

  v_scr(2,:) = dvr_mat(2,1) * v(1,:)     &
                            +            &
               dvr_mat(2,2) * v(2,:)     &
                            +            &
               dvr_mat(2,3) * v(3,:)     &
                            +            &
               dvr_mat(2,4) * v(4,:)     &
                            +            &
               dvr_mat(2,5) * v(5,:)     &
                            +            &
               dvr_mat(2,6) * v(6,:)     &
                            +            &
               dvr_mat(2,7) * v(7,:)     &
                            +            &
               dvr_mat(2,8) * v(8,:)     &
                            +            &
               dvr_mat(2,9) * v(9,:)     

  v_scr(3,:) = dvr_mat(3,1) * v(1,:)     &
                            +            &
               dvr_mat(3,2) * v(2,:)     &
                            +            &
               dvr_mat(3,3) * v(3,:)     &
                            +            &
               dvr_mat(3,4) * v(4,:)     &
                            +            &
               dvr_mat(3,5) * v(5,:)     &
                            +            &
               dvr_mat(3,6) * v(6,:)     &
                            +            &
               dvr_mat(3,7) * v(7,:)     &
                            +            &
               dvr_mat(3,8) * v(8,:)     &
                            +            &
               dvr_mat(3,9) * v(9,:)     

  v_scr(4,:) = dvr_mat(4,1) * v(1,:)     &
                            +            &
               dvr_mat(4,2) * v(2,:)     &
                            +            &
               dvr_mat(4,3) * v(3,:)     &
                            +            &
               dvr_mat(4,4) * v(4,:)     &
                            +            &
               dvr_mat(4,5) * v(5,:)     &
                            +            &
               dvr_mat(4,6) * v(6,:)     &
                            +            &
               dvr_mat(4,7) * v(7,:)     &
                            +            &
               dvr_mat(4,8) * v(8,:)     &
                            +            &
               dvr_mat(4,9) * v(9,:)     

  v_scr(5,:) = dvr_mat(5,1) * v(1,:)     &
                            +            &
               dvr_mat(5,2) * v(2,:)     &
                            +            &
               dvr_mat(5,3) * v(3,:)     &
                            +            &
               dvr_mat(5,4) * v(4,:)     &
                            +            &
               dvr_mat(5,5) * v(5,:)     &
                            +            &
               dvr_mat(5,6) * v(6,:)     &
                            +            &
               dvr_mat(5,7) * v(7,:)     &
                            +            &
               dvr_mat(5,8) * v(8,:)     &
                            +            &
               dvr_mat(5,9) * v(9,:)     

  v_scr(6,:) = dvr_mat(6,1) * v(1,:)     &
                            +            &
               dvr_mat(6,2) * v(2,:)     &
                            +            &
               dvr_mat(6,3) * v(3,:)     &
                            +            &
               dvr_mat(6,4) * v(4,:)     &
                            +            &
               dvr_mat(6,5) * v(5,:)     &
                            +            &
               dvr_mat(6,6) * v(6,:)     &
                            +            &
               dvr_mat(6,7) * v(7,:)     &
                            +            &
               dvr_mat(6,8) * v(8,:)     &
                            +            &
               dvr_mat(6,9) * v(9,:)     

  v_scr(7,:) = dvr_mat(7,1) * v(1,:)     &
                            +            &
               dvr_mat(7,2) * v(2,:)     &
                            +            &
               dvr_mat(7,3) * v(3,:)     &
                            +            &
               dvr_mat(7,4) * v(4,:)     &
                            +            &
               dvr_mat(7,5) * v(5,:)     &
                            +            &
               dvr_mat(7,6) * v(6,:)     &
                            +            &
               dvr_mat(7,7) * v(7,:)     &
                            +            &
               dvr_mat(7,8) * v(8,:)     &
                            +            &
               dvr_mat(7,9) * v(9,:)     

  v_scr(8,:) = dvr_mat(8,1) * v(1,:)     &
                            +            &
               dvr_mat(8,2) * v(2,:)     &
                            +            &
               dvr_mat(8,3) * v(3,:)     &
                            +            &
               dvr_mat(8,4) * v(4,:)     &
                            +            &
               dvr_mat(8,5) * v(5,:)     &
                            +            &
               dvr_mat(8,6) * v(6,:)     &
                            +            &
               dvr_mat(8,7) * v(7,:)     &
                            +            &
               dvr_mat(8,8) * v(8,:)     &
                            +            &
               dvr_mat(8,9) * v(9,:)     

  v_scr(9,:) = dvr_mat(9,1) * v(1,:)     &
                            +            &
               dvr_mat(9,2) * v(2,:)     &
                            +            &
               dvr_mat(9,3) * v(3,:)     &
                            +            &
               dvr_mat(9,4) * v(4,:)     &
                            +            &
               dvr_mat(9,5) * v(5,:)     &
                            +            &
               dvr_mat(9,6) * v(6,:)     &
                            +            &
               dvr_mat(9,7) * v(7,:)     &
                            +            &
               dvr_mat(9,8) * v(8,:)     &
                            +            &
               dvr_mat(9,9) * v(9,:)     
!
  v(1:9,:) = v_scr(1:9,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_9_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_10_z
!***begin prologue     v_so_mat_v_10_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_10_z
!
  SUBROUTINE v_so_mat_v_10_z(v,               &
                             v_scr,           &
                             dvr_mat,         &
                             ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)         :: v
  COMPLEX*16, DIMENSION(ni,nj)         :: v_scr
  COMPLEX*16, DIMENSION(10,10)         :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)     &
                            +            &
               dvr_mat(1,2) * v(2,:)     &
                            +            &
               dvr_mat(1,3) * v(3,:)     &
                            +            &
               dvr_mat(1,4) * v(4,:)     &
                            +            &
               dvr_mat(1,5) * v(5,:)     &
                            +            &
               dvr_mat(1,6) * v(6,:)     &
                            +            &
               dvr_mat(1,7) * v(7,:)     &
                            +            &
               dvr_mat(1,8) * v(8,:)     &
                            +            &
               dvr_mat(1,9) * v(9,:)     &
                            +            &
               dvr_mat(1,10) * v(10,:)   

  v_scr(2,:) = dvr_mat(2,1) * v(1,:)     &
                            +            &
               dvr_mat(2,2) * v(2,:)     &
                            +            &
               dvr_mat(2,3) * v(3,:)     &
                            +            &
               dvr_mat(2,4) * v(4,:)     &
                            +            &
               dvr_mat(2,5) * v(5,:)     &
                            +            &
               dvr_mat(2,6) * v(6,:)     &
                            +            &
               dvr_mat(2,7) * v(7,:)     &
                            +            &
               dvr_mat(2,8) * v(8,:)     &
                            +            &
               dvr_mat(2,9) * v(9,:)     &
                            +            &
               dvr_mat(2,10) * v(10,:)

  v_scr(3,:) = dvr_mat(3,1) * v(1,:)     &
                            +            &
               dvr_mat(3,2) * v(2,:)     &
                            +            &
               dvr_mat(3,3) * v(3,:)     &
                            +            &
               dvr_mat(3,4) * v(4,:)     &
                            +            &
               dvr_mat(3,5) * v(5,:)     &
                            +            &
               dvr_mat(3,6) * v(6,:)     &
                            +            &
               dvr_mat(3,7) * v(7,:)     &
                            +            &
               dvr_mat(3,8) * v(8,:)     &
                            +            &
               dvr_mat(3,9) * v(9,:)     &
                            +            &
               dvr_mat(3,10) * v(10,:)   

  v_scr(4,:) = dvr_mat(4,1) * v(1,:)     &
                            +            &
               dvr_mat(4,2) * v(2,:)     &
                            +            &
               dvr_mat(4,3) * v(3,:)     &
                            +            &
               dvr_mat(4,4) * v(4,:)     &
                            +            &
               dvr_mat(4,5) * v(5,:)     &
                            +            &
               dvr_mat(4,6) * v(6,:)     &
                            +            &
               dvr_mat(4,7) * v(7,:)     &
                            +            &
               dvr_mat(4,8) * v(8,:)     &
                            +            &
               dvr_mat(4,9) * v(9,:)     &
                            +            &
               dvr_mat(4,10) * v(10,:)   

  v_scr(5,:) = dvr_mat(5,1) * v(1,:)     &
                            +            &
               dvr_mat(5,2) * v(2,:)     &
                            +            &
               dvr_mat(5,3) * v(3,:)     &
                            +            &
               dvr_mat(5,4) * v(4,:)     &
                            +            &
               dvr_mat(5,5) * v(5,:)     &
                            +            &
               dvr_mat(5,6) * v(6,:)     &
                            +            &
               dvr_mat(5,7) * v(7,:)     &
                            +            &
               dvr_mat(5,8) * v(8,:)     &
                            +            &
               dvr_mat(5,9) * v(9,:)     &
                            +            &
               dvr_mat(5,10) * v(10,:)   

  v_scr(6,:) = dvr_mat(6,1) * v(1,:)     &
                            +            &
               dvr_mat(6,2) * v(2,:)     &
                            +            &
               dvr_mat(6,3) * v(3,:)     &
                            +            &
               dvr_mat(6,4) * v(4,:)     &
                            +            &
               dvr_mat(6,5) * v(5,:)     &
                            +            &
               dvr_mat(6,6) * v(6,:)     &
                            +            &
               dvr_mat(6,7) * v(7,:)     &
                            +            &
               dvr_mat(6,8) * v(8,:)     &
                            +            &
               dvr_mat(6,9) * v(9,:)     &
                            +            &
               dvr_mat(6,10) * v(10,:)   

  v_scr(7,:) = dvr_mat(7,1) * v(1,:)     &
                            +            &
               dvr_mat(7,2) * v(2,:)     &
                            +            &
               dvr_mat(7,3) * v(3,:)     &
                            +            &
               dvr_mat(7,4) * v(4,:)     &
                            +            &
               dvr_mat(7,5) * v(5,:)     &
                            +            &
               dvr_mat(7,6) * v(6,:)     &
                            +            &
               dvr_mat(7,7) * v(7,:)     &
                            +            &
               dvr_mat(7,8) * v(8,:)     &
                            +            &
               dvr_mat(7,9) * v(9,:)     &
                            +            &
               dvr_mat(7,10) * v(10,:)   

  v_scr(8,:) = dvr_mat(8,1) * v(1,:)     &
                            +            &
               dvr_mat(8,2) * v(2,:)     &
                            +            &
               dvr_mat(8,3) * v(3,:)     &
                            +            &
               dvr_mat(8,4) * v(4,:)     &
                            +            &
               dvr_mat(8,5) * v(5,:)     &
                            +            &
               dvr_mat(8,6) * v(6,:)     &
                            +            &
               dvr_mat(8,7) * v(7,:)     &
                            +            &
               dvr_mat(8,8) * v(8,:)     &
                            +            &
               dvr_mat(8,9) * v(9,:)     &
                            +            &
               dvr_mat(8,10) * v(10,:)   

  v_scr(9,:) = dvr_mat(9,1) * v(1,:)     &
                            +            &
               dvr_mat(9,2) * v(2,:)     &
                            +            &
               dvr_mat(9,3) * v(3,:)     &
                            +            &
               dvr_mat(9,4) * v(4,:)     &
                            +            &
               dvr_mat(9,5) * v(5,:)     &
                            +            &
               dvr_mat(9,6) * v(6,:)     &
                            +            &
               dvr_mat(9,7) * v(7,:)     &
                            +            &
               dvr_mat(9,8) * v(8,:)     &
                            +            &
               dvr_mat(9,9) * v(9,:)     &
                            +            &
               dvr_mat(9,10) * v(10,:)   

  v_scr(10,:) = dvr_mat(10,1) * v(1,:)   &
                              +          &
                dvr_mat(10,2) * v(2,:)   &
                              +          &
                dvr_mat(10,3) * v(3,:)   &
                              +          &
                dvr_mat(10,4) * v(4,:)   &
                              +          &
                dvr_mat(10,5) * v(5,:)   &
                              +          &
                dvr_mat(10,6) * v(6,:)   &
                              +          &
                dvr_mat(10,7) * v(7,:)   &
                              +          &
                dvr_mat(10,8) * v(8,:)   &
                              +          &
                dvr_mat(10,9) * v(9,:)   &
                              +          &
                dvr_mat(10,10) * v(10,:)   
!
  v(1:10,:) = v_scr(1:10,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_10_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_11_z
!***begin prologue     v_so_mat_v_11_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_11_z
!
  SUBROUTINE v_so_mat_v_11_z(v,          &
                             v_scr,      &
                             dvr_mat,    &
                             ni,nj)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(11,11)            :: dvr_mat
  INTEGER                                  :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) =    dvr_mat(1,1) * v(1,:)   &
                                 +        &
                  dvr_mat(1,2) * v(2,:)   &
                                 +        &
                  dvr_mat(1,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(1,4) * v(4,:)   &
                                 +        &
                  dvr_mat(1,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(1,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(1,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(1,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(1,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(1,10) * v(10,:) &   
                                 +        &
                  dvr_mat(1,11) * v(11,:)    

  v_scr(2,:) =    dvr_mat(2,1) * v(1,:)   &
                                 +        &
                  dvr_mat(2,2) * v(2,:)   &
                                 +        &
                  dvr_mat(2,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(2,4) * v(4,:)   &
                                 +        &
                  dvr_mat(2,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(2,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(2,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(2,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(2,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(2,10) * v(10,:) &   
                                 +        &
                  dvr_mat(2,11) * v(11,:)    

  v_scr(3,:) =    dvr_mat(3,1) * v(1,:)   &
                                 +        &
                  dvr_mat(3,2) * v(2,:)   &
                                 +        &
                  dvr_mat(3,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(3,4) * v(4,:)   &
                                 +        &
                  dvr_mat(3,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(3,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(3,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(3,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(3,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(3,10) * v(10,:) &   
                                 +        &
                  dvr_mat(3,11) * v(11,:)    

  v_scr(4,:) =    dvr_mat(4,1) * v(1,:)   &
                                 +        &
                  dvr_mat(4,2) * v(2,:)   &
                                 +        &
                  dvr_mat(4,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(4,4) * v(4,:)   &
                                 +        &
                  dvr_mat(4,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(4,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(4,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(4,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(4,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(4,10) * v(10,:) &   
                                 +        &
                  dvr_mat(4,11) * v(11,:)    

  v_scr(5,:) =    dvr_mat(5,1) * v(1,:)   &
                                 +        &
                  dvr_mat(5,2) * v(2,:)  &
                                 +        &
                  dvr_mat(5,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(5,4) * v(4,:)   &
                                 +        &
                  dvr_mat(5,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(5,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(5,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(5,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(5,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(5,10) * v(10,:) &   
                                 +        &
                  dvr_mat(5,11) * v(11,:)    

  v_scr(6,:) =    dvr_mat(6,1) * v(1,:)   &
                                 +        &
                  dvr_mat(6,2) * v(2,:)   &
                                 +        &
                  dvr_mat(6,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(6,4) * v(4,:)   &
                                 +        &
                  dvr_mat(6,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(6,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(6,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(6,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(6,9) * v(9,:)   &  
                                 +        &
                  dvr_mat(6,10) * v(10,:) &     
                                 +        &
                  dvr_mat(6,11) * v(11,:)    

  v_scr(7,:) =    dvr_mat(7,1) * v(1,:)   &
                                 +        &
                  dvr_mat(7,2) * v(2,:)   &
                                 +        &
                  dvr_mat(7,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(7,4) * v(4,:)   &
                                 +        &
                  dvr_mat(7,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(7,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(7,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(7,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(7,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(7,10) * v(10,:) &  
                                 +        &
                  dvr_mat(7,11) * v(11,:)    

  v_scr(8,:) =    dvr_mat(8,1) * v(1,:)   &
                                 +        &
                  dvr_mat(8,2) * v(2,:)   &
                                 +        &
                  dvr_mat(8,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(8,4) * v(4,:)   &
                                 +        &
                  dvr_mat(8,5) * v(5,:)   &  
                                 +        &
                  dvr_mat(8,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(8,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(8,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(8,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(8,10) * v(10,:) &     
                                 +        &
                  dvr_mat(9,11) * v(11,:)    

  v_scr(9,:) =    dvr_mat(9,1) * v(1,:)   &
                                 +        &
                  dvr_mat(9,2) * v(2,:)   &
                                 +        &
                  dvr_mat(9,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(9,4) * v(4,:)   &
                                 +        &
                  dvr_mat(9,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(9,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(9,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(9,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(9,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(9,10) * v(10,:) &   
                                 +        &
                  dvr_mat(9,11) * v(11,:)    

  v_scr(10,:) =   dvr_mat(10,1) * v(1,:) &
                                   +     &
                  dvr_mat(10,2) * v(2,:) &
                                   +     &
                  dvr_mat(10,3) * v(3,:) &
                                   +     &   
                  dvr_mat(10,4) * v(4,:) &
                                   +     &
                  dvr_mat(10,5) * v(5,:) &   
                                   +     &
                  dvr_mat(10,6) * v(6,:) &   
                                   +     &
                  dvr_mat(10,7) * v(7,:) &   
                                   +     &
                  dvr_mat(10,8) * v(8,:) &   
                                   +     &
                  dvr_mat(10,9) * v(9,:) &   
                                   +     &
                  dvr_mat(10,10) * v(10,:)&   
                                 +        &
                  dvr_mat(10,11) * v(11,:)     

  v_scr(11,:) =   dvr_mat(11,1) * v(1,:) &
                                   +     &
                  dvr_mat(11,2) * v(2,:) &
                                   +     &
                  dvr_mat(11,3) * v(3,:) &
                                   +     &   
                  dvr_mat(11,4) * v(4,:) &
                                   +     &
                  dvr_mat(11,5) * v(5,:) &   
                                   +     &
                  dvr_mat(11,6) * v(6,:) &   
                                   +     &
                  dvr_mat(11,7) * v(7,:) &   
                                   +     &
                  dvr_mat(11,8) * v(8,:) &   
                                   +     &
                  dvr_mat(11,9) * v(9,:) &   
                                   +     &
                  dvr_mat(11,10) * v(10,:)&   
                                 +        &
                  dvr_mat(11,11) * v(11,:)     
!
  v(1:11,:) = v_scr(1:11,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_11_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**begin prologue     v_out_v_in_so_mat
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       v_out_v_in_so_mat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck v_v_so_mat_gen_d
!***begin prologue     v_v_so_mat_gen_d
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multipiles
!***                   for the gen_deral FEDVR HAmiltonian. 
!***                   The routine is not needed in 1D problems. 
!***                   In 2D, the i index runs over the number
!                      of points in the y coordinate while in 3D, it runs over
!                      the product of the number of points in the z and y coordinates.
!
!***references
!***routines called    ebcxx, ambcxx, apbcxx
!***end prologue   
!
  SUBROUTINE v_v_so_mat_gen_d(v,          &
                              v_scr,      &
                              dvr_mat,    & 
                              ni,nj,nk)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj, nk
  INTEGER                              :: j, k
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(nk,nk)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
!
!
  v_scr(:,1:nk)= 0.d0
  DO k=1,nk
     DO j=1,nk
        v_scr(:,j) = v_scr(:,j) + v(:,k) * dvr_mat(k,j)
     END DO
  END DO
!
  v(:,1:nk) = v_scr(:,1:nk)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
END SUBROUTINE v_v_so_mat_gen_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_v_so_mat_2_d
!***begin prologue     v_v_so_mat_2_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagation matrix multiplies
!**                    for the special case of 2*2 matrices.
!
!***references
!***routines called
!***end prologue       v_v_so_mat_2_d
  SUBROUTINE v_v_so_mat_2_d(v,             &
                            v_scr,         &
                            dvr_mat,       & 
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(2,2)               :: dvr_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  v_scr(:,1) = v(:,1)  * dvr_mat(1,1)    &
                       +                 &
               v(:,2)  * dvr_mat(2,1)    
 
  v_scr(:,2) = v(:,1)  * dvr_mat(1,2)    &
                       +                 &
               v(:,2)  * dvr_mat(2,2)    

!
  v(:,1:2) = v_scr(:,1:2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_so_mat_2_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_so_mat_3_d
!***begin prologue     v_v_so_mat_3_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_v_so_mat_3_d
!
  SUBROUTINE v_v_so_mat_3_d(v,          &
                            v_scr,      &
                            dvr_mat,    & 
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(3,3)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       
  v_scr(:,1) = v(:,1)  * dvr_mat(1,1)    &
                       +                 &
               v(:,2)  * dvr_mat(2,1) &
                       +                 &
               v(:,3)  * dvr_mat(3,1)    

  v_scr(:,2) = v(:,1)  * dvr_mat(1,2)    &
                       +                 &
               v(:,2)  * dvr_mat(2,2) &
                       +                 &
               v(:,3)  * dvr_mat(3,2)    

  v_scr(:,3) = v(:,1)  * dvr_mat(1,3)    &
                       +                 &
               v(:,2)  * dvr_mat(2,3) &
                       +                 &
               v(:,3)  * dvr_mat(3,3)    

!
  v(:,1:3) = v_scr(:,1:3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_so_mat_3_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_so_mat_4_d
!***begin prologue     v_v_so_mat_4_d  
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_v_so_mat_4_d
!
  SUBROUTINE v_v_so_mat_4_d(v,        &
                            v_scr,    &
                            dvr_mat,  & 
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(4,4)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
  v_scr(:,1) = v(:,1)  * dvr_mat(1,1)    &
                       +                 &
               v(:,2)  * dvr_mat(2,1)    &
                       +                 &
               v(:,3)  * dvr_mat(3,1)    &
                       +                 &
               v(:,4)  * dvr_mat(4,1)    

  v_scr(:,2) = v(:,1)  * dvr_mat(1,2)    &
                       +                 &
               v(:,2)  * dvr_mat(2,2)    &
                       +                 &
               v(:,3)  * dvr_mat(3,2)    &
                       +                 &
               v(:,4)  * dvr_mat(4,2)    

  v_scr(:,3) = v(:,1)  * dvr_mat(1,3)    &
                       +                 &
               v(:,2)  * dvr_mat(2,3)    &
                       +                 &
               v(:,3)  * dvr_mat(3,3)    &
                       +                 &
               v(:,4)  * dvr_mat(4,3)    

  v_scr(:,4) = v(:,1)  * dvr_mat(1,4)    &
                       +                 &
               v(:,2)  * dvr_mat(2,4)    &
                       +                 &
               v(:,3)  * dvr_mat(3,4)    &
                       +                 &
               v(:,4)  * dvr_mat(4,4)    

!
  v(:,1:4) = v_scr(:,1:4)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_so_mat_4_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_so_mat_5_d
!***begin prologue     v_v_so_mat_5_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_so_mat_5_d
!
  SUBROUTINE v_v_so_mat_5_d(v,        &
                            v_scr,    &
                            dvr_mat,  & 
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(5,5)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(:,1) = v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)    

  v_scr(:,2) = v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)    

  v_scr(:,3) = v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)    

  v_scr(:,4) = v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)    

  v_scr(:,5) = v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)    

!
  v(:,1:5) = v_scr(:,1:5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_so_mat_5_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_so_mat_6_d
!***begin prologue     v_v_so_mat_6_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_so_mat_6_d
!
  SUBROUTINE v_v_so_mat_6_d(v,        &
                            v_scr,    &
                            dvr_mat,  & 
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(6,6)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
  v_scr(:,1) = v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)    

  v_scr(:,2) = v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)    

  v_scr(:,3) = v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)    

  v_scr(:,4) = v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)    

  v_scr(:,5) = v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)    

  v_scr(:,6) = v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)    

!
  v(:,1:6) = v_scr(:,1:6)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_so_mat_6_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_so_mat_7_d
!***begin prologue     v_v_so_mat_7_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_so_mat_7_d
!
  SUBROUTINE v_v_so_mat_7_d(v,        &
                            v_scr,    &
                            dvr_mat,  & 
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(7,7)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  v_scr(:,1) = v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)    

  v_scr(:,2) = v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)    

  v_scr(:,3) = v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)    

  v_scr(:,4) = v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)    

  v_scr(:,5) = v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)    

  v_scr(:,6) = v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)    

  v_scr(:,7) = v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)    

!
  v(:,1:7) = v_scr(:,1:7)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_so_mat_7_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_so_mat_8_d
!***begin prologue     v_v_so_mat_8_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_so_mat_8_d
!
  SUBROUTINE v_v_so_mat_8_d(v,        &
                            v_scr,    &
                            dvr_mat,  & 
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(8,8)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  v_scr(:,1) = v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)  &
                       +               &
               v(:,8)  * dvr_mat(8,1)    

  v_scr(:,2) = v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)    

  v_scr(:,3) = v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)    

  v_scr(:,4) = v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)  &
                       +               &
               v(:,8)  * dvr_mat(8,4)    

  v_scr(:,5) = v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)  &
                       +               &
               v(:,8)  * dvr_mat(8,5)    

  v_scr(:,6) = v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)  &
                       +               &
               v(:,8)  * dvr_mat(8,6)    

  v_scr(:,7) = v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)  &
                       +               &
               v(:,8)  * dvr_mat(8,7)    

  v_scr(:,8) = v(:,1)  * dvr_mat(1,8)  &
                       +               &
               v(:,2)  * dvr_mat(2,8)  &
                       +               &
               v(:,3)  * dvr_mat(3,8)  &
                       +               &
               v(:,4)  * dvr_mat(4,8)  &
                       +               &
               v(:,5)  * dvr_mat(5,8)  &
                       +               &
               v(:,6)  * dvr_mat(6,8)  &
                       +               &
               v(:,7)  * dvr_mat(7,8)  &
                       +               &
               v(:,8)  * dvr_mat(8,8)    

!
  v(:,1:8) = v_scr(:,1:8)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_so_mat_8_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_so_mat_9_d
!***begin prologue     v_v_so_mat_9_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_so_mat_9_d
!
  SUBROUTINE v_v_so_mat_9_d(v,        &
                            v_scr,    &
                            dvr_mat,  & 
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(9,9)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  v_scr(:,1) = v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)  &
                       +               &
               v(:,8)  * dvr_mat(8,1)  &
                       +               &
               v(:,9)  * dvr_mat(9,1)    

  v_scr(:,2) = v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)  &
                       +               &
               v(:,9)  * dvr_mat(9,2)    

  v_scr(:,3) = v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)  &
                       +               &
               v(:,8)  * dvr_mat(8,3)  &
                       +               &
               v(:,9)  * dvr_mat(9,3)    

  v_scr(:,4) = v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)  &
                       +               &
               v(:,8)  * dvr_mat(8,4)  &
                       +               &
               v(:,9)  * dvr_mat(9,4)    

  v_scr(:,5) = v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)  &
                       +               &
               v(:,8)  * dvr_mat(8,5)  &
                       +               &
               v(:,9)  * dvr_mat(9,5)    

  v_scr(:,6) = v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)  &
                       +               &
               v(:,8)  * dvr_mat(8,6)  &
                       +               &
               v(:,9)  * dvr_mat(9,6)    

  v_scr(:,7) = v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)  &
                       +               &
               v(:,8)  * dvr_mat(8,7)  &
                       +               &
               v(:,9)  * dvr_mat(9,7)    

  v_scr(:,8) = v(:,1)  * dvr_mat(1,8)  &
                       +               &
               v(:,2)  * dvr_mat(2,8)  &
                       +               &
               v(:,3)  * dvr_mat(3,8)  &
                       +               &
               v(:,4)  * dvr_mat(4,8)  &
                       +               &
               v(:,5)  * dvr_mat(5,8)  &
                       +               &
               v(:,6)  * dvr_mat(6,8)  &
                       +               &
               v(:,7)  * dvr_mat(7,8)  &
                       +               &
               v(:,8)  * dvr_mat(8,8)  &
                       +               &
               v(:,9)  * dvr_mat(9,8)    

  v_scr(:,9) = v(:,1)  * dvr_mat(1,9)  &
                       +               &
               v(:,2)  * dvr_mat(2,9)  &
                       +               &
               v(:,3)  * dvr_mat(3,9)  &
                       +               &
               v(:,4)  * dvr_mat(4,9)  &
                       +               &
               v(:,5)  * dvr_mat(5,9)  &
                       +               &
               v(:,6)  * dvr_mat(6,9)  &
                       +               &
               v(:,7)  * dvr_mat(7,9)  &
                       +               &
               v(:,8)  * dvr_mat(8,9)  &
                       +               &
               v(:,9)  * dvr_mat(9,9)    
!
  v(:,1:9) = v_scr(:,1:9)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_so_mat_9_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_so_mat_10_d
!***begin prologue     v_v_so_mat_10_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_so_mat_10_d
!
  SUBROUTINE v_v_so_mat_10_d(v,        &
                             v_scr,    &
                             dvr_mat,  & 
                             ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(ni,nj)             :: v
  REAL*8, DIMENSION(ni,nj)             :: v_scr
  REAL*8, DIMENSION(10,10)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(:,1) = v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)  &
                       +               &
               v(:,8)  * dvr_mat(8,1)  &
                       +               &
               v(:,9)  * dvr_mat(9,1)  &
                       +               &
               v(:,10) * dvr_mat(10,1)   

  v_scr(:,2) = v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)  &
                       +               &
               v(:,9)  * dvr_mat(9,2)  &
                       +               &
               v(:,10) * dvr_mat(10,2)   

  v_scr(:,3) = v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)  &
                       +               &
               v(:,8)  * dvr_mat(8,3)  &
                       +               &
               v(:,9)  * dvr_mat(9,3)  &
                       +               &
               v(:,10) * dvr_mat(10,3)   

  v_scr(:,4) = v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)  &
                       +               &
               v(:,8)  * dvr_mat(8,4)  &
                       +               &
               v(:,9)  * dvr_mat(9,4)  &
                       +               &
               v(:,10) * dvr_mat(10,4)   

  v_scr(:,5) = v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)  &
                       +               &
               v(:,8)  * dvr_mat(8,5)  &
                       +               &
               v(:,9)  * dvr_mat(9,5)  &
                       +               &
               v(:,10) * dvr_mat(10,5)   

  v_scr(:,6) = v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)  &
                       +               &
               v(:,8)  * dvr_mat(8,6)  &
                       +               &
               v(:,9)  * dvr_mat(9,6)  &
                       +               &
               v(:,10) * dvr_mat(10,6)   

  v_scr(:,7) = v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)  &
                       +               &
               v(:,8)  * dvr_mat(8,7)  &
                       +               &
               v(:,9)  * dvr_mat(9,7)  &
                       +               &
               v(:,10) * dvr_mat(10,7)   

  v_scr(:,8) = v(:,1)  * dvr_mat(1,8)  &
                       +               &
               v(:,2)  * dvr_mat(2,8)  &
                       +               &
               v(:,3)  * dvr_mat(3,8)  &
                       +               &
               v(:,4)  * dvr_mat(4,8)  &
                       +               &
               v(:,5)  * dvr_mat(5,8)  &
                       +               &
               v(:,6)  * dvr_mat(6,8)  &
                       +               &
               v(:,7)  * dvr_mat(7,8)  &
                       +               &
               v(:,8)  * dvr_mat(8,8)  &
                       +               &
               v(:,9)  * dvr_mat(9,8)  &
                       +               &
               v(:,10) * dvr_mat(10,8)   

  v_scr(:,9) = v(:,1)  * dvr_mat(1,9)  &
                       +               &
               v(:,2)  * dvr_mat(2,9)  &
                       +               &
               v(:,3)  * dvr_mat(3,9)  &
                       +               &
               v(:,4)  * dvr_mat(4,9)  &
                       +               &
               v(:,5)  * dvr_mat(5,9)  &
                       +               &
               v(:,6)  * dvr_mat(6,9)  &
                       +               &
               v(:,7)  * dvr_mat(7,9)  &
                       +               &
               v(:,8)  * dvr_mat(8,9)  &
                       +               &
               v(:,9)  * dvr_mat(9,9)  &
                       +               &
               v(:,10) * dvr_mat(10,9)   

  v_scr(:,10) = v(:,1)  * dvr_mat(1,10)  &
                       +                &
               v(:,2)  * dvr_mat(2,10)  &
                       +                &
               v(:,3)  * dvr_mat(3,10)  &
                       +                &
               v(:,4)  * dvr_mat(4,10)  &
                       +                &
               v(:,5)  * dvr_mat(5,10)  &
                       +                &
               v(:,6)  * dvr_mat(6,10)  &
                       +                &
               v(:,7)  * dvr_mat(7,10)  &
                       +                &
               v(:,8)  * dvr_mat(8,10)  &
                       +                &
               v(:,9)  * dvr_mat(9,10)  &
                       +                &
               v(:,10) * dvr_mat(10,10)   
!
  v(:,1:10) = v_scr(:,1:10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_so_mat_10_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck v_v_so_mat_gen_z
!***begin prologue     v_v_so_mat_gen_z
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multipiles
!***                   for the gen_deral FEDVR HAmiltonian. 
!***                   The routine is not needed in 1D problems. 
!***                   In 2D, the i index runs over the number
!                      of points in the y coordinate while in 3D, it runs over
!                      the product of the number of points in the z and y coordinates.
!
!***references
!***routines called    ebcxx, ambcxx, apbcxx
!***end prologue   
!
  SUBROUTINE v_v_so_mat_gen_z(v,          &
                              v_scr,      &
                              dvr_mat,    & 
                              ni,nj,nk)
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, nk
  INTEGER                                  :: j, k
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(nk,nk)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
!
!
  v_scr(:,1:nk)= (0.d0,0.d0)
  DO k=1,nk
     DO j=1,nk
        v_scr(:,j) = v_scr(:,j) + v(:,k) * dvr_mat(k,j)
     END DO
  END DO
!
  v(:,1:nk) = v_scr(:,1:nk)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
END SUBROUTINE v_v_so_mat_gen_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_so_mat_2_z
!***begin prologue     v_v_so_mat_2_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagation matrix multiplies
!**                    for the special case of 2*2 matrices.
!
!***references
!***routines called
!***end prologue       v_v_so_mat_2_z
  SUBROUTINE v_v_so_mat_2_z(v,             &
                            v_scr,         &
                            dvr_mat,       & 
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(2,2)               :: dvr_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  v_scr(:,1) = v(:,1)  * dvr_mat(1,1)    &
                       +                 &
               v(:,2)  * dvr_mat(2,1)    
 
  v_scr(:,2) = v(:,1)  * dvr_mat(1,2)    &
                       +                 &
               v(:,2)  * dvr_mat(2,2)    

!
  v(:,1:2) = v_scr(:,1:2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_so_mat_2_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_so_mat_3_z
!***begin prologue     v_v_so_mat_3_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_v_so_mat_3_z
!
  SUBROUTINE v_v_so_mat_3_z(v,          &
                            v_scr,      &
                            dvr_mat,    & 
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(3,3)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       
  v_scr(:,1) = v(:,1)  * dvr_mat(1,1)    &
                       +                 &
               v(:,2)  * dvr_mat(2,1) &
                       +                 &
               v(:,3)  * dvr_mat(3,1)    

  v_scr(:,2) = v(:,1)  * dvr_mat(1,2)    &
                       +                 &
               v(:,2)  * dvr_mat(2,2) &
                       +                 &
               v(:,3)  * dvr_mat(3,2)    

  v_scr(:,3) = v(:,1)  * dvr_mat(1,3)    &
                       +                 &
               v(:,2)  * dvr_mat(2,3) &
                       +                 &
               v(:,3)  * dvr_mat(3,3)    

!
  v(:,1:3) = v_scr(:,1:3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_so_mat_3_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_so_mat_4_z
!***begin prologue     v_v_so_mat_4_z  
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_v_so_mat_4_z
!
  SUBROUTINE v_v_so_mat_4_z(v,        &
                            v_scr,    &
                            dvr_mat,  & 
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(4,4)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
  v_scr(:,1) = v(:,1)  * dvr_mat(1,1)    &
                       +                 &
               v(:,2)  * dvr_mat(2,1)    &
                       +                 &
               v(:,3)  * dvr_mat(3,1)    &
                       +                 &
               v(:,4)  * dvr_mat(4,1)    

  v_scr(:,2) = v(:,1)  * dvr_mat(1,2)    &
                       +                 &
               v(:,2)  * dvr_mat(2,2)    &
                       +                 &
               v(:,3)  * dvr_mat(3,2)    &
                       +                 &
               v(:,4)  * dvr_mat(4,2)    

  v_scr(:,3) = v(:,1)  * dvr_mat(1,3)    &
                       +                 &
               v(:,2)  * dvr_mat(2,3)    &
                       +                 &
               v(:,3)  * dvr_mat(3,3)    &
                       +                 &
               v(:,4)  * dvr_mat(4,3)    

  v_scr(:,4) = v(:,1)  * dvr_mat(1,4)    &
                       +                 &
               v(:,2)  * dvr_mat(2,4)    &
                       +                 &
               v(:,3)  * dvr_mat(3,4)    &
                       +                 &
               v(:,4)  * dvr_mat(4,4)    

!
  v(:,1:4) = v_scr(:,1:4)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_so_mat_4_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_so_mat_5_z
!***begin prologue     v_v_so_mat_5_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_so_mat_5_z
!
  SUBROUTINE v_v_so_mat_5_z(v,        &
                            v_scr,    &
                            dvr_mat,  & 
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(5,5)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(:,1) = v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)    

  v_scr(:,2) = v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)    

  v_scr(:,3) = v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)    

  v_scr(:,4) = v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)    

  v_scr(:,5) = v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)    

!
  v(:,1:5) = v_scr(:,1:5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_so_mat_5_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_so_mat_6_z
!***begin prologue     v_v_so_mat_6_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_so_mat_6_z
!
  SUBROUTINE v_v_so_mat_6_z(v,        &
                            v_scr,    &
                            dvr_mat,  & 
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(6,6)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
  v_scr(:,1) = v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)    

  v_scr(:,2) = v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)    

  v_scr(:,3) = v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)    

  v_scr(:,4) = v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)    

  v_scr(:,5) = v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)    

  v_scr(:,6) = v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)    

!
  v(:,1:6) = v_scr(:,1:6)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_so_mat_6_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_so_mat_7_z
!***begin prologue     v_v_so_mat_7_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_so_mat_7_z
!
  SUBROUTINE v_v_so_mat_7_z(v,        &
                            v_scr,    &
                            dvr_mat,  & 
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(7,7)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  v_scr(:,1) = v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)    

  v_scr(:,2) = v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)    

  v_scr(:,3) = v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)    

  v_scr(:,4) = v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)    

  v_scr(:,5) = v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)    

  v_scr(:,6) = v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)    

  v_scr(:,7) = v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)    

!
  v(:,1:7) = v_scr(:,1:7)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_so_mat_7_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_so_mat_8_z
!***begin prologue     v_v_so_mat_8_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_so_mat_8_z
!
  SUBROUTINE v_v_so_mat_8_z(v,        &
                            v_scr,    &
                            dvr_mat,  & 
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(8,8)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  v_scr(:,1) = v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)  &
                       +               &
               v(:,8)  * dvr_mat(8,1)    

  v_scr(:,2) = v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)    

  v_scr(:,3) = v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)    

  v_scr(:,4) = v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)  &
                       +               &
               v(:,8)  * dvr_mat(8,4)    

  v_scr(:,5) = v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)  &
                       +               &
               v(:,8)  * dvr_mat(8,5)    

  v_scr(:,6) = v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)  &
                       +               &
               v(:,8)  * dvr_mat(8,6)    

  v_scr(:,7) = v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)  &
                       +               &
               v(:,8)  * dvr_mat(8,7)    

  v_scr(:,8) = v(:,1)  * dvr_mat(1,8)  &
                       +               &
               v(:,2)  * dvr_mat(2,8)  &
                       +               &
               v(:,3)  * dvr_mat(3,8)  &
                       +               &
               v(:,4)  * dvr_mat(4,8)  &
                       +               &
               v(:,5)  * dvr_mat(5,8)  &
                       +               &
               v(:,6)  * dvr_mat(6,8)  &
                       +               &
               v(:,7)  * dvr_mat(7,8)  &
                       +               &
               v(:,8)  * dvr_mat(8,8)    

!
  v(:,1:8) = v_scr(:,1:8)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_so_mat_8_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_so_mat_9_z
!***begin prologue     v_v_so_mat_9_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_so_mat_9_z
!
  SUBROUTINE v_v_so_mat_9_z(v,        &
                            v_scr,    &
                            dvr_mat,  & 
                            ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(9,9)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  v_scr(:,1) = v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)  &
                       +               &
               v(:,8)  * dvr_mat(8,1)  &
                       +               &
               v(:,9)  * dvr_mat(9,1)    

  v_scr(:,2) = v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)  &
                       +               &
               v(:,9)  * dvr_mat(9,2)    

  v_scr(:,3) = v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)  &
                       +               &
               v(:,8)  * dvr_mat(8,3)  &
                       +               &
               v(:,9)  * dvr_mat(9,3)    

  v_scr(:,4) = v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)  &
                       +               &
               v(:,8)  * dvr_mat(8,4)  &
                       +               &
               v(:,9)  * dvr_mat(9,4)    

  v_scr(:,5) = v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)  &
                       +               &
               v(:,8)  * dvr_mat(8,5)  &
                       +               &
               v(:,9)  * dvr_mat(9,5)    

  v_scr(:,6) = v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)  &
                       +               &
               v(:,8)  * dvr_mat(8,6)  &
                       +               &
               v(:,9)  * dvr_mat(9,6)    

  v_scr(:,7) = v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)  &
                       +               &
               v(:,8)  * dvr_mat(8,7)  &
                       +               &
               v(:,9)  * dvr_mat(9,7)    

  v_scr(:,8) = v(:,1)  * dvr_mat(1,8)  &
                       +               &
               v(:,2)  * dvr_mat(2,8)  &
                       +               &
               v(:,3)  * dvr_mat(3,8)  &
                       +               &
               v(:,4)  * dvr_mat(4,8)  &
                       +               &
               v(:,5)  * dvr_mat(5,8)  &
                       +               &
               v(:,6)  * dvr_mat(6,8)  &
                       +               &
               v(:,7)  * dvr_mat(7,8)  &
                       +               &
               v(:,8)  * dvr_mat(8,8)  &
                       +               &
               v(:,9)  * dvr_mat(9,8)    

  v_scr(:,9) = v(:,1)  * dvr_mat(1,9)  &
                       +               &
               v(:,2)  * dvr_mat(2,9)  &
                       +               &
               v(:,3)  * dvr_mat(3,9)  &
                       +               &
               v(:,4)  * dvr_mat(4,9)  &
                       +               &
               v(:,5)  * dvr_mat(5,9)  &
                       +               &
               v(:,6)  * dvr_mat(6,9)  &
                       +               &
               v(:,7)  * dvr_mat(7,9)  &
                       +               &
               v(:,8)  * dvr_mat(8,9)  &
                       +               &
               v(:,9)  * dvr_mat(9,9)    
!
  v(:,1:9) = v_scr(:,1:9)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_so_mat_9_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_so_mat_10_z
!***begin prologue     v_v_so_mat_10_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_so_mat_10_z
!
  SUBROUTINE v_v_so_mat_10_z(v,        &
                             v_scr,    &
                             dvr_mat,  & 
                             ni,nj)
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  COMPLEX*16, DIMENSION(10,10)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(:,1) = v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)  &
                       +               &
               v(:,8)  * dvr_mat(8,1)  &
                       +               &
               v(:,9)  * dvr_mat(9,1)  &
                       +               &
               v(:,10) * dvr_mat(10,1)   

  v_scr(:,2) = v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)  &
                       +               &
               v(:,9)  * dvr_mat(9,2)  &
                       +               &
               v(:,10) * dvr_mat(10,2)   

  v_scr(:,3) = v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)  &
                       +               &
               v(:,8)  * dvr_mat(8,3)  &
                       +               &
               v(:,9)  * dvr_mat(9,3)  &
                       +               &
               v(:,10) * dvr_mat(10,3)   

  v_scr(:,4) = v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)  &
                       +               &
               v(:,8)  * dvr_mat(8,4)  &
                       +               &
               v(:,9)  * dvr_mat(9,4)  &
                       +               &
               v(:,10) * dvr_mat(10,4)   

  v_scr(:,5) = v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)  &
                       +               &
               v(:,8)  * dvr_mat(8,5)  &
                       +               &
               v(:,9)  * dvr_mat(9,5)  &
                       +               &
               v(:,10) * dvr_mat(10,5)   

  v_scr(:,6) = v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)  &
                       +               &
               v(:,8)  * dvr_mat(8,6)  &
                       +               &
               v(:,9)  * dvr_mat(9,6)  &
                       +               &
               v(:,10) * dvr_mat(10,6)   

  v_scr(:,7) = v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)  &
                       +               &
               v(:,8)  * dvr_mat(8,7)  &
                       +               &
               v(:,9)  * dvr_mat(9,7)  &
                       +               &
               v(:,10) * dvr_mat(10,7)   

  v_scr(:,8) = v(:,1)  * dvr_mat(1,8)  &
                       +               &
               v(:,2)  * dvr_mat(2,8)  &
                       +               &
               v(:,3)  * dvr_mat(3,8)  &
                       +               &
               v(:,4)  * dvr_mat(4,8)  &
                       +               &
               v(:,5)  * dvr_mat(5,8)  &
                       +               &
               v(:,6)  * dvr_mat(6,8)  &
                       +               &
               v(:,7)  * dvr_mat(7,8)  &
                       +               &
               v(:,8)  * dvr_mat(8,8)  &
                       +               &
               v(:,9)  * dvr_mat(9,8)  &
                       +               &
               v(:,10) * dvr_mat(10,8)   

  v_scr(:,9) = v(:,1)  * dvr_mat(1,9)  &
                       +               &
               v(:,2)  * dvr_mat(2,9)  &
                       +               &
               v(:,3)  * dvr_mat(3,9)  &
                       +               &
               v(:,4)  * dvr_mat(4,9)  &
                       +               &
               v(:,5)  * dvr_mat(5,9)  &
                       +               &
               v(:,6)  * dvr_mat(6,9)  &
                       +               &
               v(:,7)  * dvr_mat(7,9)  &
                       +               &
               v(:,8)  * dvr_mat(8,9)  &
                       +               &
               v(:,9)  * dvr_mat(9,9)  &
                       +               &
               v(:,10) * dvr_mat(10,9)   

  v_scr(:,10) = v(:,1)  * dvr_mat(1,10)  &
                       +                &
               v(:,2)  * dvr_mat(2,10)  &
                       +                &
               v(:,3)  * dvr_mat(3,10)  &
                       +                &
               v(:,4)  * dvr_mat(4,10)  &
                       +                &
               v(:,5)  * dvr_mat(5,10)  &
                       +                &
               v(:,6)  * dvr_mat(6,10)  &
                       +                &
               v(:,7)  * dvr_mat(7,10)  &
                       +                &
               v(:,8)  * dvr_mat(8,10)  &
                       +                &
               v(:,9)  * dvr_mat(9,10)  &
                       +                &
               v(:,10) * dvr_mat(10,10)   
!
  v(:,1:10) = v_scr(:,1:10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_so_mat_10_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     Routines to form a packed
!                     representation of the Hamiltonian.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck pack_h.f
!***begin prologue     pack_h
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           pack, hamiltonian
!***author             schneider, b. i.(nsf)
!***source
!***purpose            The one-body matrix multiply of the DVR Hamiltonian
!                      on a vector may be performed using a packed form of
!                      the DVR hamiltonian in which only the non-zeros are
!                      used.  An arrays of indices and non-zero elements      
!                      are created in this subroutine.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       pack_h
  SUBROUTINE pack_h(dim)
  IMPLICIT NONE
  INTEGER                                :: dim
  REAL*8                                 :: eps= 1.d-15
  INTEGER                                :: posble, i, j, count
  posble=nphy(dim)*(nphy(dim)+1)/2
!
!  The diagonals are stored separately.
!
  DO  i=1,nphy(dim)
      buf(dim)%d(i) = grid(dim)%ke(i,i) 
  END DO
!
! Do the off-diagonals
! 
  count = 0
  DO  i=1,nphy(dim)
      DO  j=1,i-1
          IF(ABS(grid(dim)%ke(i,j)) > eps) THEN
             count = count + 1
             buf(dim)%hbuf(count) = grid(dim)%ke(i,j)
             buf(dim)%hibuf(1,count)=i
             buf(dim)%hibuf(2,count)=j
          END IF
      END DO
  END DO
  nonz(dim)=count
  count = count + nphy(dim)
  WRITE(iout,1) dim, count, posble
1 FORMAT (/,1X,'spatial dimension                           = ',i10,/,1x,   &
               'number of non-zero matrix elements          = ',i10,/,1x,   &
               'possible number of non-zero matrix elements = ',i10)
END SUBROUTINE pack_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck rd_hpak.f
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           pack, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            read and print a packed DVR hamiltonian
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       rd_hpak
  SUBROUTINE rd_hpak(diag,hbuf,ihbuf,n,prn)
  IMPLICIT NONE
  INTEGER                                :: n, nel, i
  REAL*8, DIMENSION(n)                   :: diag
  REAL*8, DIMENSION(*)                   :: hbuf
  INTEGER, DIMENSION(2,*)                :: ihbuf
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  CALL iosys('read real "diagonals of dvr hamiltonian" from rwf',    &
              n,diag,0,' ')
  CALL iosys('read integer "number dvr matrix elements" from rwf',1, &
              nel,0,' ')
  CALL iosys('read integer "packed dvr buffer" from rwf',2*nel,      &
              ihbuf,0,' ')
  CALL iosys('read real "packed dvr hamiltonian" from rwf',nel,      &
              hbuf,0,' ')
  IF(prn) THEN
     title='diagonals'
     CALL prntrm(title,diag,n,1,n,1,iout)
     WRITE(iout,1)
     DO  i=1,nel
         WRITE(iout,2) ihbuf(1,i), ihbuf(2,i), hbuf(i)
     END DO
  END IF
1    FORMAT(/,1X,'non-zero matrix elements',/,5X,'I',5X,'J', 5X,'Element')
2    FORMAT(i5,5X,i5,5X,e15.8)
END SUBROUTINE rd_hpak
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                     Routines to multiply a packed Hamiltonian
!                     matrix on a vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck dvr_m_v_1d_d
!***begin prologue     dvr_m_v_1d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Multiply a packed DVR hamiltonian on a vector.
!***                   This for one-dimensional, real vectors.
!***                   There are other routines which do this operation
!**                    avoiding packing and recognizing the intrisic
!**                    structure of the FEDVR Hamiltonian.  This module
!**                    was written much earlier and is based on the
!**                    typical quantum chemistry procedure of separating
!**                    out the diagonal and then packing the non-zero
!**                    elements and their indices into arrays.  The
!**                    2D array hibuf contains the indices and the array
!**                    hbuf the matrix element.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       dvr_m_v_1d_d
  SUBROUTINE dvr_m_v_1d_d(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                :: nv
  REAL*8, DIMENSION(nphy(1),nv)          :: v_in          
  REAL*8, DIMENSION(nphy(1),nv)          :: v_out          
  REAL*8                                 :: hij
  INTEGER                                :: i, j, nel 
!
!    First add in the diagonal contribution.
!
  DO i=1,nphy(1)
     v_out(i,1:nv) = v_out(i,1:nv)                       &
                                     +                   &
                     buf(1)%d(i) * v_in(i,1:nv)
  END DO
!
! Now loop over the non-zero, off-diagonal eleemnts.
!
  DO  nel=1,nonz(1)
      i=buf(1)%hibuf(1,nel)
      j=buf(1)%hibuf(2,nel)
      hij=buf(1)%hbuf(nel)
      v_out(i,1:nv) = v_out(i,1:nv) + hij * v_in(j,1:nv)
      v_out(j,1:nv) = v_out(j,1:nv) + hij * v_in(i,1:nv)
  END DO
END SUBROUTINE dvr_m_v_1d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck dvr_m_v_1d_z
!***begin prologue     dvr_m_v_1d_z
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Complex vector version of previous routine.
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       dvr_m_v_1d_z
  SUBROUTINE dvr_m_v_1d_z(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                :: nv
  COMPLEX*16, DIMENSION(nphy(1),nv)      :: v_in          
  COMPLEX*16, DIMENSION(nphy(1),nv)      :: v_out          
  REAL*8                                 :: hij
  INTEGER                                :: i, j, nel 
!
!    First add in the diagonal contribution.
!
  DO i=1,nphy(1)
     v_out(i,1:nv) = v_out(i,1:nv)                     &
                                    +                  &
                     buf(1)%d(i) * v_in(i,1:nv)
  END DO
!
! Now loop over the non-zero, off-diagonal eleemnts.
!
  DO  nel=1,nonz(1)
      i=buf(1)%hibuf(1,nel)
      j=buf(1)%hibuf(2,nel)
      hij=buf(1)%hbuf(nel)
      v_out(i,1:nv) = v_out(i,1:nv) + hij * v_in(j,1:nv)
      v_out(j,1:nv) = v_out(j,1:nv) + hij * v_in(i,1:nv)
  END DO
END SUBROUTINE dvr_m_v_1d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck dvr_m_v_2d_d.f
!***begin prologue     dvr_m_v_2d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Multiply a packed DVR hamiltonian on a vector.  This
!***                   is for two dimensional, real vectors.
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       dvr_m_v_2d_d
  SUBROUTINE dvr_m_v_2d_d(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                     :: nv
  REAL*8, DIMENSION(nphy(2),nphy(1),nv)       :: v_in, v_out
  REAL*8                                      :: hij
  INTEGER                                     :: i, nel, j
!
!                      Perform the multiply in nphy(2)
!
  DO  i=1,nphy(2)
      v_out(i,1:nphy(1),1:nv) = v_out(i,1:nphy(1),1:nv)              &
                                         +                           &
                                buf(2)%d(i) * v_in(i,1:nphy(1),1:nv)
  END DO
  DO  nel=1,nonz(2)
      i=buf(2)%hibuf(1,nel)
      j=buf(2)%hibuf(2,nel)
      hij=buf(2)%hbuf(nel)
      v_out(i,1:nphy(1),1:nv) = v_out(i,1:nphy(1),1:nv)              &
                                         +                           &
                               hij * v_in(j,1:nphy(1),1:nv)
      v_out(j,1:nphy(1),1:nv) = v_out(j,1:nphy(1),1:nv)              &
                                         +                           &
                                hij * v_in(i,1:nphy(1),1:nv)
  END DO
!
!                      Perform the multiply in nphy(1)
!
  DO  i=1,nphy(1)
      v_out(1:nphy(2),i,1:nv) = v_out(1:nphy(2),i,1:nv)              &
                                         +                           &
                                buf(1)%d(i) * v_in(1:nphy(2),i,1:nv)
  END DO
  DO  nel=1,nonz(1)
      i=buf(1)%hibuf(1,nel)
      j=buf(1)%hibuf(2,nel)
      hij=buf(1)%hbuf(nel)
      v_out(1:nphy(2),i,1:nv) = v_out(1:nphy(2),i,1:nv)              &
                                         +                           &
                               hij * v_in(1:nphy(2),j,1:nv)
      v_out(1:nphy(2),j,1:nv) = v_out(1:nphy(2),j,1:nv)              &
                                         +                           &
                                hij * v_in(1:nphy(2),i,1:nv)
  END DO
END SUBROUTINE dvr_m_v_2d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck dvr_m_v_2d_z.f
!***begin prologue     dvr_m_v_2d_z
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Same as previous for complex vectors.
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       dvr_m_v_2d_z
  SUBROUTINE dvr_m_v_2d_z(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                     :: nv
  COMPLEX*16, DIMENSION(nphy(2),nphy(1),nv)   :: v_in
  COMPLEX*16, DIMENSION(nphy(2),nphy(1),nv)   :: v_out
  REAL*8                                      :: hij
  INTEGER                                     :: i, nel, j
  DO  i=1,nphy(2)
      v_out(i,1:nphy(1),1:nv) = v_out(i,1:nphy(1),1:nv)                 &
                                         +                              &
                                buf(2)%d(i) * v_in(i,1:nphy(1),1:nv)
  END DO
  DO  nel=1,nonz(2)
      i=buf(2)%hibuf(1,nel)
      j=buf(2)%hibuf(2,nel)
      hij=buf(2)%hbuf(nel)
      v_out(i,1:nphy(1),1:nv) = v_out(i,1:nphy(1),1:nv)                  &
                                         +                               &
                                hij * v_in(j,1:nphy(1),1:nv)
      v_out(j,1:nphy(1),1:nv) = v_out(j,1:nphy(1),1:nv)                  &
                                         +                               &
                                hij * v_in(i,1:nphy(1),1:nv)
  END DO
  DO  i=1,nphy(1)
      v_out(1:nphy(2),i,1:nv) = v_out(1:nphy(2),i,1:nv)                  &
                                         +                               &
                                buf(1)%d(i) * v_in(1:nphy(2),i,1:nv)
  END DO
  DO  nel=1,nonz(1)
      i=buf(1)%hibuf(1,nel)
      j=buf(1)%hibuf(2,nel)
      hij=buf(1)%hbuf(nel)
      v_out(1:nphy(2),i,1:nv) = v_out(1:nphy(2),i,1:nv)                  &
                                               +                         &
                                hij * v_in(1:nphy(2),j,1:nv)
      v_out(1:nphy(2),j,1:nv) = v_out(1:nphy(2),j,1:nv)                  &
                                               +                         &
                                hij * v_in(1:nphy(2),i,1:nv)
  END DO
END SUBROUTINE dvr_m_v_2d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck dvr_m_v_3d_d.f
!***begin prologue     dvr_m_v_3d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Multiply a packed DVR hamiltonian on a vector.  This
!***                   is for three dimensional, real vectors.
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       dvr_m_v_3d_d
  SUBROUTINE dvr_m_v_3d_d(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                             :: nv 
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),nv)       :: v_in, v_out
  REAL*8                                              :: hij
  INTEGER                                             :: i, nel, j
!
!                Perform the multiply in nphy(3)
!
  DO i=1,nphy(3)
     v_out(i,1:nphy(2),1:nphy(1),1:nv)                     & 
                         =                                 &
     v_out(i,1:nphy(2),1:nphy(1),1:nv)                     &
                         +                                 &
     buf(3)%d(i) * v_in(i,1:nphy(2),1:nphy(1),1:nv)     
  END DO
  DO nel=1,nonz(3)
      i=buf(3)%hibuf(1,nel)
      j=buf(3)%hibuf(2,nel)
      hij=buf(3)%hbuf(nel)
      v_out(i,1:nphy(2),1:nphy(1),1:nv)                    &
                        =                                  &
      v_out(i,1:nphy(2),1:nphy(1),1:nv)                    &
                        +                                  &
           hij * v_in(j,1:nphy(2),1:nphy(1),1:nv)
      v_out(j,1:nphy(2),1:nphy(1),1:nv)                    &
                        =                                  &
      v_out(j,1:nphy(2),1:nphy(1),1:nv)                    &
                        +                                  &
      hij * v_in(i,1:nphy(2),1:nphy(1),1:nv)
  END DO
!
!                Perform the multiply in nphy(2)
!
  DO  i=1,nphy(2)
      v_out(1:nphy(3),i,1:nphy(1),1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),i,1:nphy(1),1:nv)                    &
                        +                                  &
      buf(2)%d(i) * v_in(1:nphy(3),i,1:nphy(1),1:nv) 
  END DO
  DO  nel=1,nonz(2)
      i=buf(2)%hibuf(1,nel)
      j=buf(2)%hibuf(2,nel)
      hij=buf(2)%hbuf(nel)
      v_out(1:nphy(3),i,1:nphy(2),1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),i,1:nphy(1),1:nv)                    &
                        +                                  & 
              hij * v_in(1:nphy(3),j,1:nphy(1),1:nv)
      v_out(1:nphy(3),j,1:nphy(1),1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),j,1:nphy(1),1:nv)                    &
                        +                                  &
      hij * v_in(1:nphy(3),i,1:nphy(1),1:nv)
  END DO
!
!                Perform the multiply in nphy(1)
!
  DO  i=1,nphy(1)
      v_out(1:nphy(3),1:nphy(2),i,1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),1:nphy(2),i,1:nv)                    &
                        +                                  &
      buf(1)%d(i) * v_in(1:nphy(3),1:nphy(2),i,1:nv)
  END DO
  DO  nel=1,nonz(1)
      i=buf(1)%hibuf(1,nel)
      j=buf(1)%hibuf(2,nel)
      hij=buf(1)%hbuf(nel)
      v_out(1:nphy(3),1:nphy(2),i,1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),1:nphy(2),i,1:nv)                    &
                        +                                  &
             hij * v_in(1:nphy(3),1:nphy(2),j,1:nv)        
      v_out(1:nphy(3),1:nphy(2),j,1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),1:nphy(2),j,1:nv)                    & 
                        +                                  &
      hij * v_in(1:nphy(3),1:nphy(2),i,1:nv)
  END DO
END SUBROUTINE dvr_m_v_3d_d 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck dvr_m_v_3d_z.f
!***begin prologue     dvr_m_v_3d_z
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Multiply a packed DVR hamiltonian on a vector.  Same
!***                   as previous for complex vectors. 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       dvr_m_v_3d_z
  SUBROUTINE dvr_m_v_3d_z(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                              :: nv
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1),nv)    :: v_in, v_out
  REAL*8                                               :: hij
  INTEGER                                              :: i, nel, j
  DO i=1,nphy(3)
     v_out(i,1:nphy(2),1:nphy(1),1:nv)                     & 
                         =                                 &
     v_out(i,1:nphy(2),1:nphy(1),1:nv)                     &
                         +                                 &
     buf(3)%d(i) * v_in(i,1:nphy(2),1:nphy(1),1:nv)     
  END DO
  DO nel=1,nonz(3)
      i=buf(3)%hibuf(1,nel)
      j=buf(3)%hibuf(2,nel)
      hij=buf(3)%hbuf(nel)
      v_out(i,1:nphy(2),1:nphy(1),1:nv)                    &
                        =                                  &
      v_out(i,1:nphy(2),1:nphy(1),1:nv)                    &
                        +                                  &
           hij * v_in(j,1:nphy(2),1:nphy(1),1:nv)
      v_out(j,1:nphy(2),1:nphy(1),1:nv)                    &
                        =                                  &
      v_out(j,1:nphy(2),1:nphy(1),1:nv)                    &
                        +                                  &
      hij * v_in(i,1:nphy(2),1:nphy(1),1:nv)
  END DO
  DO  i=1,nphy(2)
      v_out(1:nphy(3),i,1:nphy(1),1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),i,1:nphy(1),1:nv)                    &
                        +                                  &
      buf(2)%d(i) * v_in(1:nphy(3),i,1:nphy(1),1:nv) 
  END DO
  DO  nel=1,nonz(2)
      i=buf(2)%hibuf(1,nel)
      j=buf(2)%hibuf(2,nel)
      hij=buf(2)%hbuf(nel)
      v_out(1:nphy(3),i,1:nphy(2),1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),i,1:nphy(1),1:nv)                    &
                        +                                  & 
      hij * v_in(1:nphy(3),j,1:nphy(1),1:nv)
      v_out(1:nphy(3),j,1:nphy(1),1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),j,1:nphy(1),1:nv)                    &
                        +                                  &
      hij * v_in(1:nphy(3),i,1:nphy(1),1:nv)
  END DO
  DO  i=1,nphy(1)
      v_out(1:nphy(3),1:nphy(2),i,1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),1:nphy(2),i,1:nv)                    &
                        +                                  &
      buf(1)%d(i) * v_in(1:nphy(3),1:nphy(2),i,1:nv)
  END DO
  DO  nel=1,nonz(1)
      i=buf(1)%hibuf(1,nel)
      j=buf(1)%hibuf(2,nel)
      hij=buf(1)%hbuf(nel)
      v_out(1:nphy(3),1:nphy(2),i,1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),1:nphy(2),i,1:nv)                    &
                        +                                  &
             hij * v_in(1:nphy(3),1:nphy(2),j,1:nv)        
      v_out(1:nphy(3),1:nphy(2),j,1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),1:nphy(2),j,1:nv)                    & 
                        +                                  &
      hij * v_in(1:nphy(3),1:nphy(2),i,1:nv)
  END DO
END SUBROUTINE dvr_m_v_3d_z 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 Routines to multiply a FD, 3, 5, or 7
!                 point formula on a vector
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**begin prologue     h_on_vector
!**date written       040902   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            module for hamiltonian times vector.
!**description        this module contains all the routines necessary to
!**                   multiply a dvr Hamiltonian or a (3,5,7) point FD
!**                   Hamiltonian on a vector
!
!                     H = H(1,1) + H(2,2) + H(3,3)
!
!     Note that vectors are stored as V(nphy(3),nphy(2),nphy(1)).  Thus 
!     matrix vector multiplies in two and three dimensions must be done 
!     carefully for nphy(2) and nphy(3).  We have( repeated indices 
!     summed over), in 3D,
!    
!     V(3,2,1) =            H(1,1)   * V(3,2,1) 
!                                    + 
!                           H(2,2)   * V(3,2,1) 
!                                    + 
!                           H(3,3)   * V(3,2,1)   
!                                    =
!                           V(3,2,1) * H(1,1) 
!                                    + 
!                           V(3,2,1) * H(2,2) 
!                                    + 
!                           H(3,3) * V(3,2,1) 
!
! Thus the first mutiply may be done as a matrix V(3*2,1) on H(1,1), 
! the second as V(3,2,1) * H(2,2) with an outer loop over index 1 and 
! the third as a simple matrix vector multiply, H(3,3) * V(3,2*1).
!
!                             and in 2D,
!
!     V(2,1)   =            H(1,1) * V(2,1) 
!                                  + 
!                           H(2,2) * V(2,1)  
!                                  =
!                           V(2,1) * H(1,1) 
!                                  + 
!                           H(2,2) * V(2,1) 
!
!
!        
!**references
!**                   
!**routines called    
!**end prologue       h_on_vector
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!               The h_v and dvr_m_v routines perform the matrix
!               vector multiplication using a packed form of the matrix.
!               The finite_element_m_v routines use the element matrices
!               directly and have no need to pack the Hamiltonian.
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck h_v_d
!**begin prologue     h_v_d
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector in one dimension
!**
!**references
!**routines called h_m_v_1d_d, h_m_v_2d_d, h_m_v_3d_d
!**end prologue       h_v_d
  SUBROUTINE h_v_d(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                        :: nv, i
  REAL*8, DIMENSION(n3d,nv)      :: v_in          
  REAL*8, DIMENSION(n3d,nv)      :: v_out          
  v_out(:,:) = 0.d0
!
  IF(spdim == 1) THEN
     CALL h_m_v_1d_d(v_in,v_out,nv)
  ELSE IF(spdim == 2) THEN
     CALL h_m_v_2d_d(v_in,v_out,nv)
  ELSE IF(spdim == 3) THEN
     CALL h_m_v_3d_d(v_in,v_out,nv)
  END IF  
!
!
! take care of the diagonal v
!
  DO i=1,n3d
     v_out(i,:) = v_out(i,:) + v_tot(i) * v_in(i,:)
  END DO
  IF(log_prp(3)) THEN
     title='hamiltonian on vectors'
     CALL prntrm(title,v_out,n3d,nv,n3d,nv,iout)
  END IF
END SUBROUTINE h_v_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck h_v_z
!**begin prologue     h_v_z
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector in one dimension
!**
!**references
!**routines called h_m_v_1d_z, h_m_v_2d_z, h_m_v_3d_z
!**end prologue       h_v_z
  SUBROUTINE h_v_z(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                :: nv, i
  COMPLEX*16, DIMENSION(n3d,nv)          :: v_in          
  COMPLEX*16, DIMENSION(n3d,nv)          :: v_out          
  v_out(:,:) = (0.d0,0.d0)
  IF(spdim == 1) THEN
     CALL h_m_v_1d_z(v_in,v_out,nv)
  ELSE IF(spdim == 2) THEN
     CALL h_m_v_2d_z(v_in,v_out,nv)
  ELSE IF(spdim == 3) THEN
     CALL h_m_v_3d_z(v_in,v_out,nv)
  END IF  
!
!
! take care of the diagonal v
!
  DO i=1,n3d
     v_out(i,:) = v_out(i,:) + v_tot(i) * v_in(i,:)
  END DO
  IF(log_prp(3)) THEN
     title='hamiltonian on vectors'
     CALL prntcm(title,v_out,n3d,nv,n3d,nv,iout)
  END IF
END SUBROUTINE h_v_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck h_m_v_1d_d
!***begin prologue     h_m_v_1d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       h_m_v_1d_d
  SUBROUTINE h_m_v_1d_d(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                    :: nv
  REAL*8, DIMENSION(nphy(1),nv)              :: v_in          
  REAL*8, DIMENSION(nphy(1),nv)              :: v_out          
  REAL*8                                     :: hij
  INTEGER                                    :: i, j, nel 
  IF(typke /= 'fd') THEN
     CALL dvr_m_v_1d_d(v_in,v_out,nv)
  ELSE
     IF(row(1) == 2) THEN
        CALL ds3_bmm_d(v_in,v_out,1,1,nv,nphy(1),1)
     ELSE IF(row(1) == 3) THEN
        CALL ds5_bmm_d(v_in,v_out,1,1,nv,nphy(1),1)
     ELSE IF(row(1) == 4) THEN
        CALL ds7_bmm_d(v_in,v_out,1,1,nv,nphy(1),1)
     ELSE
        CALL lnkerr('error in band size')
     END IF
  END IF
END SUBROUTINE h_m_v_1d_d
!***********************************************************************
!deck h_m_v_1d_z
!***begin prologue     h_m_v_1d_z
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       h_m_v_1d_z
  SUBROUTINE h_m_v_1d_z(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                        :: nv
  COMPLEX*16, DIMENSION(nphy(1),nv)              :: v_in          
  COMPLEX*16, DIMENSION(nphy(1),nv)              :: v_out          
  REAL*8                                         :: hij
  INTEGER                                        :: i, j, nel 
  IF(typke /= 'fd') THEN
     CALL dvr_m_v_1d_z(v_in,v_out,nv)
  ELSE
     IF(row(1) == 2) THEN
        CALL ds3_bmm_z(v_in,v_out,1,1,nv,nphy(1),1)
     ELSE IF(row(1) == 3) THEN
        CALL ds5_bmm_z(v_in,v_out,1,1,nv,nphy(1),1)
     ELSE IF(row(1) == 4) THEN
        CALL ds7_bmm_z(v_in,v_out,1,1,nv,nphy(1),1)
     ELSE
        CALL lnkerr('error in band size')
     END IF
  END IF
END SUBROUTINE h_m_v_1d_z
!***********************************************************************
!deck h_m_v_2d_d
!***begin prologue     h_m_v_2d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       h_m_v_2d_d
  SUBROUTINE h_m_v_2d_d(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                    :: nv
  REAL*8, DIMENSION(nphy(2),nphy(1),nv)      :: v_in          
  REAL*8, DIMENSION(nphy(2),nphy(1),nv)      :: v_out          
  REAL*8                                     :: hij
  INTEGER                                    :: i, j, nel 
  IF(typke /= 'fd') THEN
     CALL dvr_m_v_2d_d(v_in,v_out,nv)
  ELSE
     IF(row(2) == 2) THEN
        CALL ds3_bmm_d(v_in,v_out,1,nv,nphy(1),nphy(2),2)
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE
           CALL lnkerr('error')
        END IF
     ELSE IF(row(2) == 3) THEN
        CALL ds5_bmm_d(v_in,v_out,1,nv,nphy(1),nphy(2),2)
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE IF(row(2) == 4) THEN
        CALL ds7_bmm_d(v_in,v_out,1,nv,nphy(1),nphy(2),2)
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE
        CALL lnkerr('error in band size')
     END IF
  END IF
END SUBROUTINE h_m_v_2d_d
!***********************************************************************
!deck h_m_v_2d_z
!***begin prologue     h_m_v_2d_z
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       h_m_v_2d_z
  SUBROUTINE h_m_v_2d_z(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                        :: nv
  COMPLEX*16, DIMENSION(nphy(2),nphy(1),nv)      :: v_in          
  COMPLEX*16, DIMENSION(nphy(2),nphy(1),nv)      :: v_out          
  REAL*8                                         :: hij
  INTEGER                                        :: i, j, nel 
  IF(typke /= 'fd') THEN
     CALL dvr_m_v_2d_z(v_in,v_out,nv)
  ELSE
     IF(row(2) == 2) THEN
        CALL ds3_bmm_z(v_in,v_out,1,nv,nphy(1),nphy(2),2)
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE
           CALL lnkerr('error')
        END IF
     ELSE IF(row(2) == 3) THEN
        CALL ds5_bmm_z(v_in,v_out,1,nv,nphy(1),nphy(2),2)
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE IF(row(2) == 4) THEN
        CALL ds7_bmm_z(v_in,v_out,1,nv,nphy(1),nphy(2),2)
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE
        CALL lnkerr('error in band size')
     END IF
  END IF
  END SUBROUTINE h_m_v_2d_z
!***********************************************************************
!deck h_m_v_3d_d
!***begin prologue     h_m_v_3d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       h_m_v_3d_d
  SUBROUTINE h_m_v_3d_d(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                            :: nv
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),nv)      :: v_in      
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),nv)      :: v_out     
  REAL*8                                             :: hij
  INTEGER                                            :: i, j, nel 
  IF(typke /= 'fd') THEN
     CALL dvr_m_v_3d_d(v_in,v_out,nv)
  ELSE
     IF(row(3) == 2) THEN
        CALL ds3_bmm_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),3)
        IF(row(2) == 2) THEN
           CALL ds3_bmmt_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 3) THEN
           CALL ds5_bmmt_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 4) THEN
           CALL ds7_bmmt_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE
           CALL lnkerr('quit')
        END IF
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE IF(row(3) == 3) THEN
        CALL ds5_bmm_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),3)
        IF(row(2) == 2) THEN
           CALL ds3_bmmt_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 3) THEN
           CALL ds5_bmmt_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 4) THEN
           CALL ds7_bmmt_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE
           CALL lnkerr('quit')
        END IF
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE IF(row(3) == 4) THEN
        CALL ds7_bmm_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),3)
        IF(row(2) == 2) THEN
           CALL ds3_bmmt_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 3) THEN
           CALL ds5_bmmt_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 4) THEN
           CALL ds7_bmmt_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE
           CALL lnkerr('quit')
        END IF
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE
        CALL lnkerr('quit')
     END IF
  END IF
END SUBROUTINE h_m_v_3d_d
!***********************************************************************
!deck h_m_v_3d_z
!***begin prologue     h_m_v_3d_z
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       h_m_v_3d_z
  SUBROUTINE h_m_v_3d_z(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                                :: nv
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1),nv)      :: v_in  
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1),nv)      :: v_out  
  REAL*8                                                 :: hij
  INTEGER                                                :: i, j 
  INTEGER                                                :: nel 
  IF(typke /= 'fd') THEN
     CALL dvr_m_v_3d_z(v_in,v_out,nv)
  ELSE
     IF(row(3) == 2) THEN
        CALL ds3_bmm_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),3)
        IF(row(2) == 2) THEN
           CALL ds3_bmmt_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 3) THEN
           CALL ds5_bmmt_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 4) THEN
           CALL ds7_bmmt_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE
           CALL lnkerr('quit')
        END IF
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE IF(row(3) == 3) THEN
        CALL ds5_bmm_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),3)
        IF(row(2) == 2) THEN
           CALL ds3_bmmt_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 3) THEN
           CALL ds5_bmmt_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 4) THEN
           CALL ds7_bmmt_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE
           CALL lnkerr('quit')
        END IF
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE IF(row(3) == 4) THEN
        CALL ds7_bmm_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),3)
        IF(row(2) == 2) THEN
           CALL ds3_bmmt_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 3) THEN
           CALL ds5_bmmt_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 4) THEN
           CALL ds7_bmmt_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE
           CALL lnkerr('quit')
        END IF
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE
        CALL lnkerr('quit')
     END IF
  END IF
  END SUBROUTINE h_m_v_3d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck ds3_bmm_d
!**begin prologue     ds3_bmm_d
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            Three point FD Hamiltonian times space vector
!**
!**description        the two non-zero elements of the symmetric
!**                   tridiagonal matrix are stored by columns as
!\begin{eqnarray}
!              m(i,i) &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds3_bmm_d
  SUBROUTINE ds3_bmm_d(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  REAL*8, DIMENSION(n4,n3,n2,n1)         :: v_in, v_out
  INTEGER                                :: i
!
!\center In the tridiagonal case all but the first and last rows
!\center "feel" all three elements of the matrix

  DO  i=2,n4-1
      
!\begin{eqnarray}
! v_{out}(i,j,k) &=&  v_{out}(i,j,k) + m(i,i+1) v_{in}(i+1,j,k)
!                                    + m(i,i)   v_{in}(i,j,k)
!                                    + m(i,i-1) v_{in}(i-1,j,k) \nonumber \\
!            &=&  v_{out}(i,j,k)     + m(i+1,i) v_{in}(i+1,j,k) \nonumber
!                                    + m(i,i)   v_{in}(i,j,k)
!                                    + m(i,i-1) v_{in}(i-1,j,k)
!\end{eqnarray}
!    Since the element (i-1) is not explictly stored in column i, we
!    use the off-diagonal element from column (i-1)
      
      v_out(i,:,:,:) = v_out(i,:,:,:)                            &
                                   +                             &
                       grid(dim)%ke(2,i) * v_in(i+1,:,:,:)       &
                                   +                             &
                       grid(dim)%ke(1,i) * v_in(i,:,:,:)         &
                                   +                             &
                       grid(dim)%ke(2,i-1) * v_in(i-1,:,:,:)
  END DO

!\center Do the special case of the first and last row.

      v_out(1,:,:,:)  = v_out(1,:,:,:)                           & 
                                   +                             &
                        grid(dim)%ke(2,1) * v_in(2,:,:,:)        &
                                   +                             &
                        grid(dim)%ke(1,1) * v_in(1,:,:,:)
      v_out(n4,:,:,:) = v_out(n4,:,:,:)                          &
                                   +                             &
                        grid(dim)%ke(2,n4-1) * v_in(n4-1,:,:,:)  & 
                                   +                             &
                        grid(dim)%ke(1,n4)   * v_in(n4,:,:,:)
END SUBROUTINE ds3_bmm_d
!***********************************************************************
!deck ds3_bmm_z
!**begin prologue     ds3_bmm_z
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            Same as previous for complex vector.
!**
!**description        the two non-zero elements of the symmetric
!**                   tridiagonal matrix are stored by columns as
!\begin{eqnarray}
!              m(i,i) &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds3_bmm_z
  SUBROUTINE ds3_bmm_z(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  COMPLEX*16, DIMENSION(n4,n3,n2,n1)     :: v_in, v_out
  INTEGER                                :: i
!
!\center In the tridiagonal case all but the first and last rows
!\center "feel" all three elements of the matrix

  DO  i=2,n4-1
      
!\begin{eqnarray}
! v_{out}(i,j,k) &=&  v_{out}(i,j,k) + m(i,i+1) v_{in}(i+1,j,k)
!                                    + m(i,i)   v_{in}(i,j,k)
!                                    + m(i,i-1) v_{in}(i-1,j,k) \nonumber \\
!            &=&  v_{out}(i,j,k)     + m(i+1,i) v_{in}(i+1,j,k) \nonumber
!                                    + m(i,i)   v_{in}(i,j,k)
!                                    + m(i,i-1) v_{in}(i-1,j,k)
!\end{eqnarray}
!    Since the element (i-1) is not explictly stored in column i, we
!    use the off-diagonal element from column (i-1)
      
      v_out(i,:,:,:) = v_out(i,:,:,:)                            &
                                   +                             &
                       grid(dim)%ke(2,i) * v_in(i+1,:,:,:)       &
                                   +                             &
                       grid(dim)%ke(1,i) * v_in(i,:,:,:)         &
                                   +                             &
                       grid(dim)%ke(2,i-1) * v_in(i-1,:,:,:)
  END DO

!\center Do the special case of the first and last row.

      v_out(1,:,:,:)  = v_out(1,:,:,:)                           & 
                                   +                             &
                        grid(dim)%ke(2,1) * v_in(2,:,:,:)        &
                                   +                             &
                        grid(dim)%ke(1,1) * v_in(1,:,:,:)
      v_out(n4,:,:,:) = v_out(n4,:,:,:)                          &
                                   +                             &
                        grid(dim)%ke(2,n4-1) * v_in(n4-1,:,:,:)  & 
                                   +                             &
                        grid(dim)%ke(1,n4)   * v_in(n4,:,:,:)
END SUBROUTINE ds3_bmm_z
!***********************************************************************
!deck ds3_bmmt_d.f
!**begin prologue     ds3_bmmt_d
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**
!**description        Same as previous but for transpose of Hamiltonian.
!\begin{eqnarray}
!              m(i,i) &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds3_bmmt_d
  SUBROUTINE ds3_bmmt_d(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  REAL*8, DIMENSION(n4,n3,n2,n1)         :: v_in, v_out
  INTEGER                                :: i
!\center In the tridiagonal case all but the first and last rows
!\center "feel" all three elements of the matrix

  DO  i=2,n3-1
      
!\begin{eqnarray}
! v_{out}(j,i,k) &=&  v_{out}(j,i,k) + m(i,i+1) v_{in}(j,i+1,k)
!                                    + m(i,i)   v_{in}(j,i,k)
!                                    + m(i,i-1) v_{in}(j,i-1,k) \nonumber \\
!              &=&  v_{out}(j,i,k)   + m(i+1,i) v_{in}(j,i+1,k) \nonumber
!                                    + m(i,i)   v_{in}(j,i,k)
!                                    + m(i,i-1) v_{in}(j,i-1,k)
!\end{eqnarray}
!    Since the element (i-1) is not explictly stored in column i, we
!    use the off-diagonal element from column (i-1)
      
      v_out(:,i,:,:) = v_out(:,i,:,:)                              &
                                         +                         &
                       grid(dim)%ke(2,i)   * v_in(:,i+1,:,:)       &
                                         +                         &
                       grid(dim)%ke(1,i)   * v_in(:,i,:,:)         &
                                         +                         &
                       grid(dim)%ke(2,i-1) * v_in(:,i-1,:,:)
  END DO

!\center Do the special case of the first and last row.

  v_out(:,1,:,:)  = v_out(:,1,:,:)                                 &
                                         +                         &
                    grid(dim)%ke(2,1)    * v_in(:,2,:,:)           &
                                         +                         &
                    grid(dim)%ke(1,1)    * v_in(:,1,:,:)
  v_out(:,n3,:,:) = v_out(:,n3,:,:)                                &
                                         +                         &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3-1,:,:)        &
                                         +                         &
                    grid(dim)%ke(1,n3)   * v_in(:,n3,:,:)
END SUBROUTINE ds3_bmmt_d
!***********************************************************************
!deck ds3_bmmt_z.f
!**begin prologue     ds3_bmmt_z
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**
!**description        See previous.
!\begin{eqnarray}
!              m(i,i) &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds3_bmmt_z
  SUBROUTINE ds3_bmmt_z(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  COMPLEX*16, DIMENSION(n4,n3,n2,n1)     :: v_in, v_out
  INTEGER                                :: i
!\center In the tridiagonal case all but the first and last rows
!\center "feel" all three elements of the matrix

  DO  i=2,n3-1
      
!\begin{eqnarray}
! v_{out}(j,i,k) &=&  v_{out}(j,i,k) + m(i,i+1) v_{in}(j,i+1,k)
!                                    + m(i,i)   v_{in}(j,i,k)
!                                    + m(i,i-1) v_{in}(j,i-1,k) \nonumber \\
!              &=&  v_{out}(j,i,k)   + m(i+1,i) v_{in}(j,i+1,k) \nonumber
!                                    + m(i,i)   v_{in}(j,i,k)
!                                    + m(i,i-1) v_{in}(j,i-1,k)
!\end{eqnarray}
!    Since the element (i-1) is not explictly stored in column i, we
!    use the off-diagonal element from column (i-1)
      
      v_out(:,i,:,:) = v_out(:,i,:,:)                              &
                                         +                         &
                       grid(dim)%ke(2,i)   * v_in(:,i+1,:,:)       &
                                         +                         &
                       grid(dim)%ke(1,i)   * v_in(:,i,:,:)         &
                                         +                         &
                       grid(dim)%ke(2,i-1) * v_in(:,i-1,:,:)
  END DO

!\center Do the special case of the first and last row.

  v_out(:,1,:,:)  = v_out(:,1,:,:)                                 &
                                         +                         &
                    grid(dim)%ke(2,1)    * v_in(:,2,:,:)           &
                                         +                         &
                    grid(dim)%ke(1,1)    * v_in(:,1,:,:)
  v_out(:,n3,:,:) = v_out(:,n3,:,:)                                &
                                         +                         &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3-1,:,:)        &
                                         +                         &
                    grid(dim)%ke(1,n3)   * v_in(:,n3,:,:)
END SUBROUTINE ds3_bmmt_z
!***********************************************************************
!deck ds5_bmm_d
!**begin prologue     ds5_bmm_d
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            Five point FD hamiltonian times space vector.
!**description        the three non-zero elements of the pentadiagonal
!**                   matrix are stored by columns as

!\begin{eqnarray}
!              m(i,i)   &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& hb(3,i) = m(i,i+2)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds5_bmm_d
  SUBROUTINE ds5_bmm_d(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                            :: n1, n2, n3, n4, dim
  REAL*8, DIMENSION(n4,n3,n2,n1)     :: v_in, v_out
  INTEGER                            :: i
!
!\center In the pentadiagonal case all but the first two and last two rows
!\center "feel" all three elements of the matrix

  DO  i=3,n4-2
!\begin{eqnarray}
! v_{out}(i,j,k) =  v_{out}(i,j,k) + m(i,i+2) v_{in}(i+2,j,k)
!                                  + m(i,i+1) v_{in}(i+1,j,k)
!                                  + m(i,i)   v_{in}(i,j,k)  \nonumber \\
!                                  + m(i,i-2) v_{in}(i-2,j,k)
!                                  + m(i,i-1) v_{in}(i-1,j,k) \nonumber \\
!                =  v_{out}(i,j,k) + m(i+2,i) v_{in}(i+2,j,k)
!                                  + m(i+1,i) v_{in}(i+1,j,k)
!                                  + m(i,i)   v_{in}(i,j,k) \nonumber \\
!                                  + m(i,i-1) v_{in}(i-1,j,k)
!                                  + m(i,i-2) v_{in}(i-2,j,k)
!\end{eqnarray}
!    Since the elements (i-1) and (i-2) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1) and (i-2)
      
      v_out(i,:,:,:) = v_out(i,:,:,:)                                    &
                                      +                                  &
                       grid(dim)%ke(3,i)   * v_in(i+2,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(2,i)   * v_in(i+1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(1,i)   * v_in(i,:,:,:)               &
                                      +                                  &
                       grid(dim)%ke(2,i-1) * v_in(i-1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(3,i-2) * v_in(i-2,:,:,:)
  END DO
!
!\center Do the special cases of the first two and last two rows.
!
  v_out(1,:,:,:)   =  v_out(1,:,:,:)                                    &
                                      +                                 &
                      grid(dim)%ke(1,1) * v_in(1,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(2,1) * v_in(2,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(3,1) * v_in(3,:,:,:)
  v_out(2,:,:,:)   =  v_out(2,:,:,:)                                    &
                                      +                                 &
                      grid(dim)%ke(2,1) * v_in(1,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(1,2) * v_in(2,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(2,2) * v_in(3,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(3,2) * v_in(4,:,:,:)
  v_out(n4-1,:,:,:) = v_out(n4-1,:,:,:)                                 &
                                      +                                 &
                      grid(dim)%ke(2,n4-1) * v_in(n4,:,:,:)             &
                                      +                                 &
                      grid(dim)%ke(1,n4-1) * v_in(n4-1,:,:,:)           &
                                      +                                 & 
                      grid(dim)%ke(2,n4-2) * v_in(n4-2,:,:,:)           &
                                      +                                 &
                      grid(dim)%ke(3,n4-3) * v_in(n4-3,:,:,:)
  v_out(n4,:,:,:)   = v_out(n4,:,:,:)                                   &
                                      +                                 &
                      grid(dim)%ke(1,n4)   * v_in(n4,:,:,:)             &
                                      +                                 &
                      grid(dim)%ke(2,n4-1) * v_in(n4-1,:,:,:)           & 
                                      +                                 &
                      grid(dim)%ke(3,n4-2) * v_in(n4-2,:,:,:)
END SUBROUTINE ds5_bmm_d
!***********************************************************************
!deck ds5_bmm_z
!**begin prologue     ds5_bmm_z
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            See previous
!**description        the three non-zero elements of the pentadiagonal
!**                   matrix are stored by columns as

!\begin{eqnarray}
!              m(i,i)   &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& hb(3,i) = m(i,i+2)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds5_bmm_z
  SUBROUTINE ds5_bmm_z(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  COMPLEX*16, DIMENSION(n4,n3,n2,n1)     :: v_in, v_out
  INTEGER                                :: i
!
!\center In the pentadiagonal case all but the first two and last two rows
!\center "feel" all three elements of the matrix

  DO  i=3,n4-2
!\begin{eqnarray}
! v_{out}(i,j,k) =  v_{out}(i,j,k) + m(i,i+2) v_{in}(i+2,j,k)
!                                  + m(i,i+1) v_{in}(i+1,j,k)
!                                  + m(i,i)   v_{in}(i,j,k)  \nonumber \\
!                                  + m(i,i-2) v_{in}(i-2,j,k)
!                                  + m(i,i-1) v_{in}(i-1,j,k) \nonumber \\
!                =  v_{out}(i,j,k) + m(i+2,i) v_{in}(i+2,j,k)
!                                  + m(i+1,i) v_{in}(i+1,j,k)
!                                  + m(i,i)   v_{in}(i,j,k) \nonumber \\
!                                  + m(i,i-1) v_{in}(i-1,j,k)
!                                  + m(i,i-2) v_{in}(i-2,j,k)
!\end{eqnarray}
!    Since the elements (i-1) and (i-2) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1) and (i-2)
      
      v_out(i,:,:,:) = v_out(i,:,:,:)                                    &
                                      +                                  &
                       grid(dim)%ke(3,i)   * v_in(i+2,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(2,i)   * v_in(i+1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(1,i)   * v_in(i,:,:,:)               &
                                      +                                  &
                       grid(dim)%ke(2,i-1) * v_in(i-1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(3,i-2) * v_in(i-2,:,:,:)
  END DO
!
!\center Do the special cases of the first two and last two rows.
!
  v_out(1,:,:,:)   =  v_out(1,:,:,:)                                    &
                                      +                                 &
                      grid(dim)%ke(1,1) * v_in(1,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(2,1) * v_in(2,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(3,1) * v_in(3,:,:,:)
  v_out(2,:,:,:)   =  v_out(2,:,:,:)                                    &
                                      +                                 &
                      grid(dim)%ke(2,1) * v_in(1,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(1,2) * v_in(2,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(2,2) * v_in(3,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(3,2) * v_in(4,:,:,:)
  v_out(n4-1,:,:,:) = v_out(n4-1,:,:,:)                                 &
                                      +                                 &
                      grid(dim)%ke(2,n4-1) * v_in(n4,:,:,:)             &
                                      +                                 &
                      grid(dim)%ke(1,n4-1) * v_in(n4-1,:,:,:)           &
                                      +                                 & 
                      grid(dim)%ke(2,n4-2) * v_in(n4-2,:,:,:)           &
                                      +                                 &
                      grid(dim)%ke(3,n4-3) * v_in(n4-3,:,:,:)
  v_out(n4,:,:,:)   = v_out(n4,:,:,:)                                   &
                                      +                                 &
                      grid(dim)%ke(1,n4)   * v_in(n4,:,:,:)             &
                                      +                                 &
                      grid(dim)%ke(2,n4-1) * v_in(n4-1,:,:,:)           & 
                                      +                                 &
                      grid(dim)%ke(3,n4-2) * v_in(n4-2,:,:,:)
END SUBROUTINE ds5_bmm_z
!***********************************************************************
!deck ds5_bmmt_d.f
!**begin prologue     ds5_bmmt_d
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            Same as previous but for transpose of hamiltonian
!**description        the three non-zero elements of the pentadiagonal
!**                   matrix are stored by columns as
!\begin{eqnarray}
!              m(i,i)   &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& hb(3,i) = m(i,i+2)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds5_bmmt_d
  SUBROUTINE ds5_bmmt_d(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  REAL*8,     DIMENSION(n4,n3,n2,n1)        :: v_in, v_out
  INTEGER                                :: i
!\center In the pentadiagonal case all but the first two and last two rows
!\center "feel" all three elements of the matrix

  DO i=3,n3-2
      
!\begin{eqnarray}
! v_{out}(j,i,k) =  v_{out}(j,i,k) + m(i,i+2) v_{in}(j,i+2,k)
!                                  + m(i,i+1) v_{in}(j,i+1,k)
!                                  + m(i,i)   v_{in}(j,i,k)  \nonumber \\
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!                                  + m(i,i-1) v_{in}(j,i-1,k) \nonumber \\
!                =  v_{out}(j,i,k) + m(i+2,i) v_{in}(j,i+2,k)
!                                  + m(i+1,i) v_{in}(j,i+1,k)
!                                  + m(i,i)   v_{in}(j,i,k) \nonumber \\
!                                  + m(i,i-1) v_{in}(j,i-1,k)
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!\end{eqnarray}
!    Since the elements (i-1) and (i-2) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1) and (i-2)
      
     v_out(:,i,:,:) = v_out(:,i,:,:)                              &  
                                  +                               &
                    grid(dim)%ke(3,i)   * v_in(:,i+2,:,:)         &
                                  +                               &
                    grid(dim)%ke(2,i)   * v_in(:,i+1,:,:)         &
                                  +                               &
                    grid(dim)%ke(1,i)   * v_in(:,i,:,:)           &
                                  +                               &
                    grid(dim)%ke(2,i-1) * v_in(:,i-1,:,:)         &
                                  +                               &
                    grid(dim)%ke(3,i-2) * v_in(:,i-2,:,:)
  END DO

!\center Do the special cases of the first two and last two rows.

  v_out(:,1,:,:)  = v_out(:,1,:,:)                                &
                                  +                               &
                  grid(dim)%ke(1,1) * v_in(:,1,:,:)               &
                                  +                               &
                  grid(dim)%ke(2,1) * v_in(:,2,:,:)               &
                                  +                               &
                  grid(dim)%ke(3,1) * v_in(:,3,:,:)
  v_out(:,2,:,:)  = v_out(:,2,:,:)                                &
                                  +                               & 
                  grid(dim)%ke(2,1) * v_in(:,1,:,:)               &
                                  +                               &
                  grid(dim)%ke(1,2) * v_in(:,2,:,:)               &
                                  +                               &
                  grid(dim)%ke(2,2) * v_in(:,3,:,:)               &
                                  +                               &
                  grid(dim)%ke(3,2) * v_in(:,4,:,:)
  v_out(:,n3-1,:,:) = v_out(:,n3-1,:,:)                           &
                                  +                               &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3,:,:)         &
                                  +                               &
                    grid(dim)%ke(1,n3-1) * v_in(:,n3-1,:,:)       &
                                  +                               &
                    grid(dim)%ke(2,n3-2) * v_in(:,n3-2,:,:)       &
                                  +                               &
                    grid(dim)%ke(3,n3-3) * v_in(:,n3-3,:,:)
  v_out(:,n3,:,:)   = v_out(:,n3,:,:)                             &
                                  +                               &
                    grid(dim)%ke(1,n3)   * v_in(:,n3,:,:)         &
                                  +                               &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3-1,:,:)       &
                                  +                               &
                    grid(dim)%ke(3,n3-2) * v_in(:,n3-2,:,:)
END SUBROUTINE ds5_bmmt_d
!***********************************************************************
!deck ds5_bmmt_z.f
!**begin prologue     ds5_bmmt_z
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            See previous
!**description        the three non-zero elements of the pentadiagonal
!**                   matrix are stored by columns as
!\begin{eqnarray}
!              m(i,i)   &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& hb(3,i) = m(i,i+2)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds5_bmmt_z
  SUBROUTINE ds5_bmmt_z(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                 :: n1, n2, n3, n4, dim
  COMPLEX*16, DIMENSION(n4,n3,n2,n1)      :: v_in, v_out
  INTEGER                                 :: i
!\center In the pentadiagonal case all but the first two and last two rows
!\center "feel" all three elements of the matrix

  DO i=3,n3-2
      
!\begin{eqnarray}
! v_{out}(j,i,k) =  v_{out}(j,i,k) + m(i,i+2) v_{in}(j,i+2,k)
!                                  + m(i,i+1) v_{in}(j,i+1,k)
!                                  + m(i,i)   v_{in}(j,i,k)  \nonumber \\
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!                                  + m(i,i-1) v_{in}(j,i-1,k) \nonumber \\
!                =  v_{out}(j,i,k) + m(i+2,i) v_{in}(j,i+2,k)
!                                  + m(i+1,i) v_{in}(j,i+1,k)
!                                  + m(i,i)   v_{in}(j,i,k) \nonumber \\
!                                  + m(i,i-1) v_{in}(j,i-1,k)
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!\end{eqnarray}
!    Since the elements (i-1) and (i-2) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1) and (i-2)
      
     v_out(:,i,:,:) = v_out(:,i,:,:)                              &  
                                  +                               &
                    grid(dim)%ke(3,i)   * v_in(:,i+2,:,:)         &
                                  +                               &
                    grid(dim)%ke(2,i)   * v_in(:,i+1,:,:)         &
                                  +                               &
                    grid(dim)%ke(1,i)   * v_in(:,i,:,:)           &
                                  +                               &
                    grid(dim)%ke(2,i-1) * v_in(:,i-1,:,:)         &
                                  +                               &
                    grid(dim)%ke(3,i-2) * v_in(:,i-2,:,:)
  END DO

!\center Do the special cases of the first two and last two rows.

  v_out(:,1,:,:)  = v_out(:,1,:,:)                                &
                                  +                               &
                  grid(dim)%ke(1,1) * v_in(:,1,:,:)               &
                                  +                               &
                  grid(dim)%ke(2,1) * v_in(:,2,:,:)               &
                                  +                               &
                  grid(dim)%ke(3,1) * v_in(:,3,:,:)
  v_out(:,2,:,:)  = v_out(:,2,:,:)                                &
                                  +                               & 
                  grid(dim)%ke(2,1) * v_in(:,1,:,:)               &
                                  +                               &
                  grid(dim)%ke(1,2) * v_in(:,2,:,:)               &
                                  +                               &
                  grid(dim)%ke(2,2) * v_in(:,3,:,:)               &
                                  +                               &
                  grid(dim)%ke(3,2) * v_in(:,4,:,:)
  v_out(:,n3-1,:,:) = v_out(:,n3-1,:,:)                           &
                                  +                               &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3,:,:)         &
                                  +                               &
                    grid(dim)%ke(1,n3-1) * v_in(:,n3-1,:,:)       &
                                  +                               &
                    grid(dim)%ke(2,n3-2) * v_in(:,n3-2,:,:)       &
                                  +                               &
                    grid(dim)%ke(3,n3-3) * v_in(:,n3-3,:,:)
  v_out(:,n3,:,:)   = v_out(:,n3,:,:)                             &
                                  +                               &
                    grid(dim)%ke(1,n3)   * v_in(:,n3,:,:)         &
                                  +                               &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3-1,:,:)       &
                                  +                               &
                    grid(dim)%ke(3,n3-2) * v_in(:,n3-2,:,:)
END SUBROUTINE ds5_bmmt_z
!***********************************************************************
!deck ds7_bmm_d.f
!**begin prologue     ds7_bmm_d
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            Seven point FD hamiltonian times space vector
!**description        the four non-zero elements of the heptadiagonal
!**                   matrix are stored by columns as

!\begin{eqnarray}
!              m(i,i) &=& grid(dim)%(1,i) \nonumber \\
!              m(i+1,i) &=& grid(dim)%(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& grid(dim)%(3,i) = m(i,i+2) \nonumber \\
!              m(i+3,i) &=& grid(dim)%(4,i) = m(i,i+3)
!\end{eqnarray}

!**references
!**routines called
!**end prologue       ds7_bmm_d
  SUBROUTINE ds7_bmm_d(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  REAL*8, DIMENSION(n4,n3,n2,n1)         :: v_in, v_out
  INTEGER                                :: i
!
!\center In the heptadiagonal case all but the first three and last three rows
!\center "feel" all four elements of the matrix
!
  DO  i=4,n4-3
      
!\begin{eqnarray}
! v_{out}(i,:,:) =  v_{out}(i,:,:) + m(i,i+3) v_{in}(i+3,:,:)
!                                  + m(i,i+2) v_{in}(i+2,:,:)
!                                  + m(i,i+1) v_{in}(i+1,:,:)  \nonumber \\
!                                  + m(i,i)   v_{in}(i,:,:)
!                                  + m(i,i-1) v_{in}(i-1,:,:)
!                                  + m(i,i-2) v_{in}(i-2,:,:)
!                                  + m(i,i-3) v_{in}(i-3,:,:) \nonumber \\
!                =  v_{out}(i,:,:) + m(i+3,i) v_{in}(i+3,:,:)
!                                  + m(i+2,i) v_{in}(i+2,:,:)
!                                  + m(i+1,i) v_{in}(i+1,:,:) \nonumber \\
!                                  + m(i,i)   v_{in}(i,:,:)
!                                  + m(i,i-1) v_{in}(i-1,:,:)
!                                  + m(i,i-2) v_{in}(i-2,:,:)
!                                  + m(i,i-3) v_{in}(i-3,:,:)
!\end{eqnarray}
!    Since the elements (i-1), (i-2) and (i-3) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1), (i-2) and (i-3)

      v_out(i,:,:,:) = v_out(i,:,:,:)                                    &
                                      +                                  &
                       grid(dim)%ke(4,i)   * v_in(i+3,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(3,i)   * v_in(i+2,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(2,i)   * v_in(i+1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(1,i)   * v_in(i,:,:,:)               &
                                      +                                  &
                       grid(dim)%ke(2,i-1) * v_in(i-1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(3,i-2) * v_in(i-2,:,:,:)             &
                                      +                                  &   
                       grid(dim)%ke(4,i-3) * v_in(i-3,:,:,:)               
  END DO

!\center Do the special cases of the first three and last three rows.

    v_out(1,:,:,:) = v_out(1,:,:,:)                                      &
                                      +                                  &
                   grid(dim)%ke(1,1) * v_in(1,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,1) * v_in(2,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(3,1) * v_in(3,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(4,1) * v_in(4,:,:,:)
    v_out(2,:,:,:) = v_out(2,:,:,:)                                      &
                                      +                                  &
                   grid(dim)%ke(2,1) * v_in(1,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(1,2) * v_in(2,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,2) * v_in(3,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(3,2) * v_in(4,:,:,:)                     & 
                                      +                                  &
                   grid(dim)%ke(4,2) * v_in(5,:,:,:)
    v_out(3,:,:,:) = v_out(3,:,:,:)                                      &
                                      +                                  &
                   grid(dim)%ke(3,1) * v_in(1,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,2) * v_in(2,:,:,:)                     & 
                                      +                                  &
                   grid(dim)%ke(1,3) * v_in(3,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,3) * v_in(4,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(3,3) * v_in(5,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(4,3) * v_in(6,:,:,:)
    v_out(n4-2,:,:,:) = v_out(n4-2,:,:,:)                                &
                                      +                                  &
                      grid(dim)%ke(3,n4-2) * v_in(n4,:,:,:)              &
                                      +                                  &
                      grid(dim)%ke(2,n4-2) * v_in(n4-1,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(1,n4-2) * v_in(n4-2,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(2,n4-3) * v_in(n4-3,:,:,:)            & 
                                      +                                  &
                      grid(dim)%ke(3,n4-4) * v_in(n4-4,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(4,n4-5) * v_in(n4-5,:,:,:)
    v_out(n4-1,:,:,:) = v_out(n4-1,:,:,:)                                &
                                      +                                  &
                      grid(dim)%ke(2,n4-1) * v_in(n4,:,:,:)              &
                                      +                                  &
                      grid(dim)%ke(1,n4-1) * v_in(n4-1,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(2,n4-2) * v_in(n4-2,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(3,n4-3) * v_in(n4-3,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(4,n4-4) * v_in(n4-4,:,:,:)
    v_out(n4,:,:,:)   = v_out(n4,:,:,:)                                  &
                                      +                                  &
                      grid(dim)%ke(1,n4)   * v_in(n4,:,:,:)              &
                                      +                                  &
                      grid(dim)%ke(2,n4-1) * v_in(n4-1,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(3,n4-2) * v_in(n4-2,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(4,n4-3) * v_in(n4-3,:,:,:)
END SUBROUTINE ds7_bmm_d
!***********************************************************************
!!deck ds7_bmm_z.f
!**begin prologue     ds7_bmm_z
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            See previous
!**description        the four non-zero elements of the heptadiagonal
!**                   matrix are stored by columns as

!\begin{eqnarray}
!              m(i,i) &=& grid(dim)%(1,i) \nonumber \\
!              m(i+1,i) &=& grid(dim)%(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& grid(dim)%(3,i) = m(i,i+2) \nonumber \\
!              m(i+3,i) &=& grid(dim)%(4,i) = m(i,i+3)
!\end{eqnarray}

!**references
!**routines called
!**end prologue       ds7_bmm_z
  SUBROUTINE ds7_bmm_z(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  COMPLEX*16, DIMENSION(n4,n3,n2,n1)     :: v_in, v_out
  INTEGER                                :: i
!
!\center In the heptadiagonal case all but the first three and last three rows
!\center "feel" all four elements of the matrix
!
  DO  i=4,n4-3
      
!\begin{eqnarray}
! v_{out}(i,:,:) =  v_{out}(i,:,:) + m(i,i+3) v_{in}(i+3,:,:)
!                                  + m(i,i+2) v_{in}(i+2,:,:)
!                                  + m(i,i+1) v_{in}(i+1,:,:)  \nonumber \\
!                                  + m(i,i)   v_{in}(i,:,:)
!                                  + m(i,i-1) v_{in}(i-1,:,:)
!                                  + m(i,i-2) v_{in}(i-2,:,:)
!                                  + m(i,i-3) v_{in}(i-3,:,:) \nonumber \\
!                =  v_{out}(i,:,:) + m(i+3,i) v_{in}(i+3,:,:)
!                                  + m(i+2,i) v_{in}(i+2,:,:)
!                                  + m(i+1,i) v_{in}(i+1,:,:) \nonumber \\
!                                  + m(i,i)   v_{in}(i,:,:)
!                                  + m(i,i-1) v_{in}(i-1,:,:)
!                                  + m(i,i-2) v_{in}(i-2,:,:)
!                                  + m(i,i-3) v_{in}(i-3,:,:)
!\end{eqnarray}
!    Since the elements (i-1), (i-2) and (i-3) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1), (i-2) and (i-3)

      v_out(i,:,:,:) = v_out(i,:,:,:)                                    &
                                      +                                  &
                       grid(dim)%ke(4,i)   * v_in(i+3,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(3,i)   * v_in(i+2,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(2,i)   * v_in(i+1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(1,i)   * v_in(i,:,:,:)               &
                                      +                                  &
                       grid(dim)%ke(2,i-1) * v_in(i-1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(3,i-2) * v_in(i-2,:,:,:)             &
                                      +                                  &   
                       grid(dim)%ke(4,i-3) * v_in(i-3,:,:,:)               
  END DO

!\center Do the special cases of the first three and last three rows.

    v_out(1,:,:,:) = v_out(1,:,:,:)                                      &
                                      +                                  &
                   grid(dim)%ke(1,1) * v_in(1,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,1) * v_in(2,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(3,1) * v_in(3,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(4,1) * v_in(4,:,:,:)
    v_out(2,:,:,:) = v_out(2,:,:,:)                                      &
                                      +                                  &
                   grid(dim)%ke(2,1) * v_in(1,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(1,2) * v_in(2,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,2) * v_in(3,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(3,2) * v_in(4,:,:,:)                     & 
                                      +                                  &
                   grid(dim)%ke(4,2) * v_in(5,:,:,:)
    v_out(3,:,:,:) = v_out(3,:,:,:)                                      &
                                      +                                  &
                   grid(dim)%ke(3,1) * v_in(1,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,2) * v_in(2,:,:,:)                     & 
                                      +                                  &
                   grid(dim)%ke(1,3) * v_in(3,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,3) * v_in(4,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(3,3) * v_in(5,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(4,3) * v_in(6,:,:,:)
    v_out(n4-2,:,:,:) = v_out(n4-2,:,:,:)                                &
                                      +                                  &
                      grid(dim)%ke(3,n4-2) * v_in(n4,:,:,:)              &
                                      +                                  &
                      grid(dim)%ke(2,n4-2) * v_in(n4-1,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(1,n4-2) * v_in(n4-2,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(2,n4-3) * v_in(n4-3,:,:,:)            & 
                                      +                                  &
                      grid(dim)%ke(3,n4-4) * v_in(n4-4,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(4,n4-5) * v_in(n4-5,:,:,:)
    v_out(n4-1,:,:,:) = v_out(n4-1,:,:,:)                                &
                                      +                                  &
                      grid(dim)%ke(2,n4-1) * v_in(n4,:,:,:)              &
                                      +                                  &
                      grid(dim)%ke(1,n4-1) * v_in(n4-1,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(2,n4-2) * v_in(n4-2,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(3,n4-3) * v_in(n4-3,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(4,n4-4) * v_in(n4-4,:,:,:)
    v_out(n4,:,:,:)   = v_out(n4,:,:,:)                                  &
                                      +                                  &
                      grid(dim)%ke(1,n4)   * v_in(n4,:,:,:)              &
                                      +                                  &
                      grid(dim)%ke(2,n4-1) * v_in(n4-1,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(3,n4-2) * v_in(n4-2,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(4,n4-3) * v_in(n4-3,:,:,:)
END SUBROUTINE ds7_bmm_z
!***********************************************************************
!deck ds7_bmmt_d.f
!**begin prologue     ds7_bmmt_d
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            Same as previous for transpose of hamiltonian
!**description        the four non-zero elements of the heptadiagonal
!**                   matrix are stored by columns as
!\begin{eqnarray}
!              m(i,i) &=& grid(dim)%ke(1,i) \nonumber \\
!              m(i+1,i) &=& grid(dim)%ke(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& grid(dim)%ke(3,i) = m(i,i+2) \nonumber \\
!              m(i+3,i) &=& grid(dim)%ke(4,i) = m(i,i+3)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds7_bmmt_d
  SUBROUTINE ds7_bmmt_d(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  REAL*8, DIMENSION(n4,n3,n2,n1)         :: v_in, v_out
  INTEGER                                :: i
!\center In the heptadiagonal case all but the first three and last three rows
!\center "feel" all four elements of the matrix

  DO i=4,n3-3
      
!\begin{eqnarray}
! v_{out}(j,i,k) =  v_{out}(j,i,k) + m(i,i+3) v_{in}(j,i+3,k)
!                                  + m(i,i+2) v_{in}(j,i+2,k)
!                                  + m(i,i+1) v_{in}(j,i+1,k) \nonumber \\
!                                  + m(i,i)   v_{in}(j,i,k)
!                                  + m(i,i-1) v_{in}(j,i-1,k)
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!                                  + m(i,i-3) v_{in}(j,i-3,k) \nonumber \\
!                =  v_{out}(j,i,k) + m(i+3,i) v_{in}(j,i+3,k)
!                                  + m(i+2,i) v_{in}(j,i+2,k)
!                                  + m(i+1,i) v_{in}(j,i+1,k) \nonumber \\
!                                  + m(i,i)   v_{in}(j,i,k)
!                                  + m(i,i-1) v_{in}(j,i-1,k)
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!                                  + m(i,i-3) v_{in}(j,i-3,k)
!\end{eqnarray}
!    Since the elements (i-1), (i-2) and (i-3) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1), (i-2) and (i-3)
      
     v_out(:,i,:,:) = v_out(:,i,:,:)                                &
                                      +                             &
                    grid(dim)%ke(4,i)   * v_in(:,i+3,:,:)           &
                                      +                             &
                    grid(dim)%ke(3,i)   * v_in(:,i+2,:,:)           & 
                                      +                             &
                    grid(dim)%ke(2,i)   * v_in(:,i+1,:,:)           &
                                      +                             &
                    grid(dim)%ke(1,i)   * v_in(:,i,:,:)             &
                                      +                             &
                    grid(dim)%ke(4,i-3) * v_in(:,i-3,:,:)           &
                                      +                             &
                    grid(dim)%ke(3,i-2) * v_in(:,i-2,:,:)           &
                                      +                             &
                    grid(dim)%ke(2,i-1) * v_in(:,i-1,:,:)
  END DO

!\center Do the special cases of the first three and last three rows.

  v_out(:,1,:,:) = v_out(:,1,:,:)                                   &
                                   +                                &
                 grid(dim)%ke(1,1) * v_in(:,1,:,:)                  &
                                   +                                &
                 grid(dim)%ke(2,1) * v_in(:,2,:,:)                  &
                                   +                                &
                 grid(dim)%ke(3,1) * v_in(:,3,:,:)                  &
                                   +                                &
                 grid(dim)%ke(4,1) * v_in(:,4,:,:)
  v_out(:,2,:,:) = v_out(:,2,:,:)                                   &
                                   +                                &
                 grid(dim)%ke(2,1) * v_in(:,1,:,:)                  & 
                                   +                                &
                 grid(dim)%ke(1,2) * v_in(:,2,:,:)                  & 
                                   +                                &
                 grid(dim)%ke(2,2) * v_in(:,3,:,:)                  &
                                   +                                &
                 grid(dim)%ke(3,2) * v_in(:,4,:,:)                  &
                                   +                                &
                 grid(dim)%ke(4,2) * v_in(:,5,:,:)
  v_out(:,3,:,:) = v_out(:,3,:,:)                                   &
                                   +                                &
                 grid(dim)%ke(3,1) * v_in(:,1,:,:)                  &
                                   +                                &
                 grid(dim)%ke(2,2) * v_in(:,2,:,:)                  &
                                   +                                &
                 grid(dim)%ke(1,3) * v_in(:,3,:,:)                  &
                                   +                                &
                 grid(dim)%ke(2,3) * v_in(:,4,:,:)                  &
                                   +                                &
                 grid(dim)%ke(3,3) * v_in(:,5,:,:)                  &
                                   +                                &
                 grid(dim)%ke(4,3) * v_in(:,6,:,:)
  v_out(:,n3-2,:,:) = v_out(:,n3-2,:,:)                             &
                                   +                                &
                    grid(dim)%ke(3,n3-2) * v_in(:,n3,:,:)           &
                                   +                                &
                   grid(dim)%ke(2,n3-2) * v_in(:,n3-1,:,:)          &
                                   +                                &
                   grid(dim)%ke(1,n3-2) * v_in(:,n3-2,:,:)          &
                                   +                                &
                   grid(dim)%ke(2,n3-3) * v_in(:,n3-3,:,:)          &
                                   +                                &
                   grid(dim)%ke(3,n3-4) * v_in(:,n3-4,:,:)          &
                                   +                                &
                   grid(dim)%ke(4,n3-5) * v_in(:,n3-5,:,:)
  v_out(:,n3-1,:,:) = v_out(:,n3-1,:,:)                             &
                                   +                                &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3,:,:)           &
                                   +                                &
                    grid(dim)%ke(1,n3-1) * v_in(:,n3-1,:,:)         &
                                   +                                &
                    grid(dim)%ke(2,n3-2) * v_in(:,n3-2,:,:)         &
                                   +                                &
                    grid(dim)%ke(3,n3-3) * v_in(:,n3-3,:,:)         &
                                   +                                &
                    grid(dim)%ke(4,n3-4) * v_in(:,n3-4,:,:)
  v_out(:,n3,:,:)   = v_out(:,n3,:,:)                               &
                                   +                                &
                    grid(dim)%ke(1,n3)   * v_in(:,n3,:,:)           &
                                   +                                &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3-1,:,:)         &
                                   +                                &
                    grid(dim)%ke(3,n3-2) * v_in(:,n3-2,:,:)         &
                                   +                                &
                    grid(dim)%ke(4,n3-3) * v_in(:,n3-3,:,:)
END SUBROUTINE ds7_bmmt_d
!***********************************************************************
!deck ds7_bmmt_z.f
!**begin prologue     ds7_bmmt_z
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            See previous
!**description        the four non-zero elements of the heptadiagonal
!**                   matrix are stored by columns as
!\begin{eqnarray}
!              m(i,i) &=& grid(dim)%ke(1,i) \nonumber \\
!              m(i+1,i) &=& grid(dim)%ke(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& grid(dim)%ke(3,i) = m(i,i+2) \nonumber \\
!              m(i+3,i) &=& grid(dim)%ke(4,i) = m(i,i+3)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds7_bmmt_z
  SUBROUTINE ds7_bmmt_z(v_in,v_out,n1,n2,n3,n4,dim)
  IMPLICIT NONE
  INTEGER                                  :: n1, n2, n3, n4, dim
  COMPLEX*16, DIMENSION(n4,n3,n2,n1)       :: v_in, v_out
  INTEGER                                  :: i
!\center In the heptadiagonal case all but the first three and last three rows
!\center "feel" all four elements of the matrix

  DO i=4,n3-3
      
!\begin{eqnarray}
! v_{out}(j,i,k) =  v_{out}(j,i,k) + m(i,i+3) v_{in}(j,i+3,k)
!                                  + m(i,i+2) v_{in}(j,i+2,k)
!                                  + m(i,i+1) v_{in}(j,i+1,k) \nonumber \\
!                                  + m(i,i)   v_{in}(j,i,k)
!                                  + m(i,i-1) v_{in}(j,i-1,k)
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!                                  + m(i,i-3) v_{in}(j,i-3,k) \nonumber \\
!                =  v_{out}(j,i,k) + m(i+3,i) v_{in}(j,i+3,k)
!                                  + m(i+2,i) v_{in}(j,i+2,k)
!                                  + m(i+1,i) v_{in}(j,i+1,k) \nonumber \\
!                                  + m(i,i)   v_{in}(j,i,k)
!                                  + m(i,i-1) v_{in}(j,i-1,k)
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!                                  + m(i,i-3) v_{in}(j,i-3,k)
!\end{eqnarray}
!    Since the elements (i-1), (i-2) and (i-3) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1), (i-2) and (i-3)
      
     v_out(:,i,:,:) = v_out(:,i,:,:)                                &
                                      +                             &
                    grid(dim)%ke(4,i)   * v_in(:,i+3,:,:)           &
                                      +                             &
                    grid(dim)%ke(3,i)   * v_in(:,i+2,:,:)           & 
                                      +                             &
                    grid(dim)%ke(2,i)   * v_in(:,i+1,:,:)           &
                                      +                             &
                    grid(dim)%ke(1,i)   * v_in(:,i,:,:)             &
                                      +                             &
                    grid(dim)%ke(4,i-3) * v_in(:,i-3,:,:)           &
                                      +                             &
                    grid(dim)%ke(3,i-2) * v_in(:,i-2,:,:)           &
                                      +                             &
                    grid(dim)%ke(2,i-1) * v_in(:,i-1,:,:)
  END DO

!\center Do the special cases of the first three and last three rows.

  v_out(:,1,:,:) = v_out(:,1,:,:)                                   &
                                   +                                &
                 grid(dim)%ke(1,1) * v_in(:,1,:,:)                  &
                                   +                                &
                 grid(dim)%ke(2,1) * v_in(:,2,:,:)                  &
                                   +                                &
                 grid(dim)%ke(3,1) * v_in(:,3,:,:)                  &
                                   +                                &
                 grid(dim)%ke(4,1) * v_in(:,4,:,:)
  v_out(:,2,:,:) = v_out(:,2,:,:)                                   &
                                   +                                &
                 grid(dim)%ke(2,1) * v_in(:,1,:,:)                  & 
                                   +                                &
                 grid(dim)%ke(1,2) * v_in(:,2,:,:)                  & 
                                   +                                &
                 grid(dim)%ke(2,2) * v_in(:,3,:,:)                  &
                                   +                                &
                 grid(dim)%ke(3,2) * v_in(:,4,:,:)                  &
                                   +                                &
                 grid(dim)%ke(4,2) * v_in(:,5,:,:)
  v_out(:,3,:,:) = v_out(:,3,:,:)                                   &
                                   +                                &
                 grid(dim)%ke(3,1) * v_in(:,1,:,:)                  &
                                   +                                &
                 grid(dim)%ke(2,2) * v_in(:,2,:,:)                  &
                                   +                                &
                 grid(dim)%ke(1,3) * v_in(:,3,:,:)                  &
                                   +                                &
                 grid(dim)%ke(2,3) * v_in(:,4,:,:)                  &
                                   +                                &
                 grid(dim)%ke(3,3) * v_in(:,5,:,:)                  &
                                   +                                &
                 grid(dim)%ke(4,3) * v_in(:,6,:,:)
  v_out(:,n3-2,:,:) = v_out(:,n3-2,:,:)                             &
                                   +                                &
                    grid(dim)%ke(3,n3-2) * v_in(:,n3,:,:)           &
                                   +                                &
                   grid(dim)%ke(2,n3-2) * v_in(:,n3-1,:,:)          &
                                   +                                &
                   grid(dim)%ke(1,n3-2) * v_in(:,n3-2,:,:)          &
                                   +                                &
                   grid(dim)%ke(2,n3-3) * v_in(:,n3-3,:,:)          &
                                   +                                &
                   grid(dim)%ke(3,n3-4) * v_in(:,n3-4,:,:)          &
                                   +                                &
                   grid(dim)%ke(4,n3-5) * v_in(:,n3-5,:,:)
  v_out(:,n3-1,:,:) = v_out(:,n3-1,:,:)                             &
                                   +                                &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3,:,:)           &
                                   +                                &
                    grid(dim)%ke(1,n3-1) * v_in(:,n3-1,:,:)         &
                                   +                                &
                    grid(dim)%ke(2,n3-2) * v_in(:,n3-2,:,:)         &
                                   +                                &
                    grid(dim)%ke(3,n3-3) * v_in(:,n3-3,:,:)         &
                                   +                                &
                    grid(dim)%ke(4,n3-4) * v_in(:,n3-4,:,:)
  v_out(:,n3,:,:)   = v_out(:,n3,:,:)                               &
                                   +                                &
                    grid(dim)%ke(1,n3)   * v_in(:,n3,:,:)           &
                                   +                                &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3-1,:,:)         &
                                   +                                &
                    grid(dim)%ke(3,n3-2) * v_in(:,n3-2,:,:)         &
                                   +                                &
                    grid(dim)%ke(4,n3-3) * v_in(:,n3-3,:,:)
END SUBROUTINE ds7_bmmt_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            END       MODULE matrix_vector_multiply_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
