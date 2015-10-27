!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    MODULE finite_element_matrix_multiply
                    INTERFACE finite_element_m_v
                MODULE PROCEDURE finite_element_m_v_d,     &
                                 finite_element_m_v_z
                    END INTERFACE finite_element_m_v
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
                          CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
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
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  USE dvr_matrix_vector_multiply
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
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  USE dvr_matrix_vector_multiply
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
END MODULE finite_element_matrix_multiply

