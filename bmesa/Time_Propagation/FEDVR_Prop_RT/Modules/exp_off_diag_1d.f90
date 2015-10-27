  SUBROUTINE exp_off_diag_1d(wave_function,scratch_vector)
     USE dvrprop_global_rt
     USE dvr_shared
     USE dvr_global
     USE exp_off_diagonal
     IMPLICIT NONE
     REAL*8, DIMENSION(nphy(1),2)                   :: wave_function
     REAL*8, DIMENSION(nphy(1),2)                   :: scratch_vector

      call exp_off_diag_m_v(wave_function(1,1),          &
                            wave_function(1,2),          &
                            scratch_vector(1,1),         &
                            scratch_vector(1,2),         &
                            nphy(1),1,1)
  END SUBROUTINE exp_off_diag_1d
