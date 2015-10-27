  SUBROUTINE exp_off_diag_3d(wave_function,scratch_vector)
     USE dvrprop_global_rt
     USE dvr_shared
     USE dvr_global
     USE exp_off_diagonal
     USE plot_wavefunction
     IMPLICIT NONE
     REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: wave_function
     REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: scratch_vector

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            This code is for the general problem involving nphy(3),
!            H(3,3) * V(3,2,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      call exp_off_diag_m_v(wave_function(1,1,1,1),         &
                            wave_function(1,1,1,2),         &
                            scratch_vector(1,1,1,1),        &
                            scratch_vector(1,1,1,2),        &
                            nphy(3),nphy(2)*nphy(1),3)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general case involving the matrix multiply 
!        involving nphy(2), that is, V(3,2,1) * H(2,2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call exp_off_diag_v_m_3_d                             &
                        (wave_function(1,1,1,1),            &
                         wave_function(1,1,1,2),            &
                         scratch_vector(1,1,1,1),           &
                         scratch_vector(1,1,1,2),           &
                         nphy(3),nphy(2),nphy(1),2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general problem involving the matrix multiply 
!        for nphy(1), that is, V(3,2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call exp_off_diag_v_m_2_d                            &
                       (wave_function(1,1,1,1),            &
                        wave_function(1,1,1,2),            &
                        scratch_vector(1,1,1,1),           &
                        scratch_vector(1,1,1,2),           &
                        nphy(3)*nphy(2),nphy(1),1)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END SUBROUTINE exp_off_diag_3d
