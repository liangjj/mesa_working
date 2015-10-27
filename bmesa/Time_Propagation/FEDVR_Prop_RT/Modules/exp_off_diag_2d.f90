  SUBROUTINE exp_off_diag_2d(wave_function,scratch_vector)
     USE dvrprop_global_rt
     USE dvr_shared
     USE dvr_global
     USE exp_off_diagonal
     IMPLICIT NONE
     REAL*8, DIMENSION(nphy(2),nphy(1),2)           :: wave_function
     REAL*8, DIMENSION(nphy(2),nphy(1),2)           :: scratch_vector
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!            This code is for the general problem involving nphy(2),
!            that is H(2,2) * V(2,1) 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      call exp_off_diag_m_v(wave_function(1,1,1),          &
                            wave_function(1,1,2),          &
                            scratch_vector(1,1,1),         &
                            scratch_vector(1,1,2),         &
                            nphy(2),nphy(1),2)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!            This code is for the general problem involving nphy(1),
!            V(2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      call exp_off_diag_v_m_2_d(wave_function(1,1,1),      &
                                wave_function(1,1,2),      &
                                scratch_vector(1,1,1),     &
                                scratch_vector(1,1,2),     &
                                nphy(2),nphy(1),1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END SUBROUTINE exp_off_diag_2d
