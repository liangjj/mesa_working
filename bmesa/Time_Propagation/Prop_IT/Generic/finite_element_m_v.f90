!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             MODULE finite_element_m_v
                         USE dvrprop_global_it
                         USE dvrprop_global_rt
                         USE dvr_shared
                         USE dvr_global
                         USE dvr_mat_mul
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck finite_element_m_v_1_d_d
!***begin prologue     finite_element_m_v_1_d_d    
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
!***end prologue       finite_element_m_v_1d_d                     
!
  SUBROUTINE finite_element_m_v_1_d_d(wave_function,scratch_vector)
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1))               :: wave_function
  REAL*8, DIMENSION(nphy(1))               :: scratch_vector
  CALL dvr_mat_mul_d(wave_function,scratch_vector,nphy(1),1,1)
  END SUBROUTINE finite_element_m_v_1_d_d
!
  SUBROUTINE finite_element_m_v_1_d_z(wave_function,scratch_vector)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(1))                :: wave_function
  COMPLEX*16, DIMENSION(nphy(1))                :: scratch_vector
  CALL dvr_mat_mul_z(wave_function,scratch_vector,nphy(1),1,1)
  END SUBROUTINE finite_element_m_v_1_d_z
!
  SUBROUTINE finite_element_m_v_2_d_d(wave_function,scratch_vector)
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1))           :: wave_function
  REAL*8, DIMENSION(nphy(2),nphy(1))           :: scratch_vector
  INTEGER                                      :: i
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       This code is for the general problem involving nphy(2),
!       that is H(2,2) * V(2,1) 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL dvr_mat_mul_d(wave_function,           &
                           scratch_vector,          &
                           nphy(2),nphy(1),2)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       This code is for the general problem involving nphy(1),
!       V(2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        CALL dvr_mat_mul_2_d(wave_function,      &
                             scratch_vector,     &
                             nphy(2),nphy(1),1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END SUBROUTINE finite_element_m_v_2_d_d
!
  SUBROUTINE finite_element_m_v_2_d_z(wave_function,scratch_vector)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(2),nphy(1))           :: wave_function
  COMPLEX*16, DIMENSION(nphy(2),nphy(1))           :: scratch_vector
  INTEGER                                           :: i
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       This code is for the general problem involving nphy(2),
!       that is H(2,2) * V(2,1) 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        CALL dvr_mat_mul_z(wave_function,scratch_vector,       &
                           nphy(2),nphy(1),2)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       This code is for the general problem involving nphy(1),
!       V(2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        CALL dvr_mat_mul_2_z(wave_function,                     &
                             scratch_vector,                    &
                             nphy(2),nphy(1),1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END SUBROUTINE finite_element_m_v_2_d_z
!
  SUBROUTINE finite_element_m_v_3_d_d(wave_function,scratch_vector)
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))        :: wave_function
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))        :: scratch_vector
  INTEGER                                           :: i, j

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            This code is for the general problem involving nphy(3),
!            H(3,3) * V(3,2,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     CALL dvr_mat_mul_d(wave_function,scratch_vector,nphy(3),     &
                        nphy(2)*nphy(1),3)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general case involving the matrix multiply 
!        involving nphy(2), that is, V(3,2,1) * H(2,2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     DO i=1,nphy(1)
        CALL dvr_mat_mul_2_d(wave_function(:,:,i),           &
                             scratch_vector(:,:,i),          &
                             nphy(3),nphy(2),2)
     END DO
   
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general problem involving the matrix multiply 
!        for nphy(1), that is, V(3,2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     CALL dvr_mat_mul_2_d(wave_function,              &
                          scratch_vector,             &
                          nphy(3)*nphy(2),nphy(1),1)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
END SUBROUTINE finite_element_m_v_3_d_d
  SUBROUTINE finite_element_m_v_3_d_z(wave_function,scratch_vector)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1))   :: wave_function
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1))   :: scratch_vector
  INTEGER                                          :: i
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            This code is for the general problem involving nphy(3),
!            H(3,3) * V(3,2,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     CALL dvr_mat_mul_z(wave_function,scratch_vector,nphy(3),     &
                        nphy(2)*nphy(1),3)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general case involving the matrix multiply 
!        involving nphy(2), that is, V(3,2,1) * H(2,2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     DO i=1,nphy(1)
        CALL dvr_mat_mul_2_z(wave_function(1,1,i),                &
                             scratch_vector(1,1,i),               &
                             nphy(3),nphy(2),2)
     END DO
   
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general problem involving the matrix multiply 
!        for nphy(1), that is, V(3,2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     CALL dvr_mat_mul_2_z(wave_function,                          &
                          scratch_vector,                         &
                          nphy(3)*nphy(2),nphy(1),1)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
END SUBROUTINE finite_element_m_v_3_d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          END MODULE finite_element_m_v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
