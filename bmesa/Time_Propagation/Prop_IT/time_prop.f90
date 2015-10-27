!deck time_prop.f
!***begin prologue     time_prop
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           split-operator, propagation, RSP2, RSP4
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate the sector time propagators
!***                   based on Lie-Trotter-Suziki formulation.
!***                   A second and fourth order method have
!***                   been programmed.  the splitting is based
!***                   on the finite element DVR approach.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       time_prop
  Subroutine time_prop
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  USE finite_element_propagator
  USE psi_h_psi
  USE h_on_vector
  IMPLICIT NONE
  REAL*4                                   :: secnds
  INTEGER                                  :: i, j
  INTEGER, DIMENSION(2)                    :: words
!
!        Calculate the time dependent propagators
!
  DO i=1,spdim
     WRITE(iout,1) i
     ALLOCATE(exp_tmp_d(maxmem_reg,maxmem_reg))
     p_fac = 1.d0
     n_prop =1
     IF(prop_order == 4) THEN
        p_fac = 1.d0/( 4.d0 - 4.d0**(1.d0/3.d0) )
        n_prop = 2
     END IF
     DO j=1,num_reg(i)
        ALLOCATE(mat_reg_d(j,i)%exp_t_mat                     &
                 (nfun_reg(j,i),nfun_reg(j,i),n_prop))
     END DO
     CALL propagator(i)
     DEALLOCATE(exp_tmp_d)
     DO j=1,num_reg(i)
        DEALLOCATE(mat_reg_d(j,i)%ke_mat_d,                   &
                    mat_reg_d(j,i)%eigvec_mat_d,              &
                    mat_reg_d(j,i)%eigval_mat_d)
     END DO
  END DO
1 FORMAT(/,15x,'Sector Time Propagtors Variable = ',i2)
END SUBROUTINE time_prop
