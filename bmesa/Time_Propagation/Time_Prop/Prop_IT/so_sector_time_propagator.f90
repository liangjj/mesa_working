!deck SO_Sector_Time_Propagator
!***begin prologue     Split_Operator_Sector_Time_Propagator
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
!***end prologue       SO_Sector_Time_Propagator
  Subroutine SO_Sector_Time_Propagator
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  USE regional_propagators
  USE h_on_vector
  IMPLICIT NONE
  REAL*4                                   :: secnds
  INTEGER                                  :: i, j
!
!        Calculate the  time dependent propagators
!
  ALLOCATE(exp_d(maxmem_reg,maxmem_reg))
  DO i=1,spdim
     WRITE(iout,1) i
     p_fac = 1.d0
     n_prop =1
     IF(prop_order == 4) THEN
        p_fac = 1.d0/( 4.d0 - 4.d0**(1.d0/3.d0) )
        n_prop = 2
     END IF
     DO j=1,num_reg(i)
        ALLOCATE(mat_reg(j,i)%exp_t_mat_d                        &
                 (nfun_reg(j,i),nfun_reg(j,i),n_prop))         
     END DO
     starting_reg = 1
     ending_reg = num_reg(i)
     n_reg = num_reg(i)
     write(iout,2)
     call reg_prop_d(i,deltat)
     DO j=1,num_reg(i)
        DEALLOCATE(mat_reg(j,i)%ke_mat_d,                      &
                   mat_reg(j,i)%eigval_mat_d,                  &
                   mat_reg(j,i)%eigvec_mat_d)
     END DO
  END DO
  DEALLOCATE(exp_d)
1 FORMAT(/,15x,'Sector Time Propagtors Variable = ',i2)
2 FORMAT('Constructing Propagators for Sectors with Real '     &
            'Potentials')
3 FORMAT('Constructing Propagators for Sectors with Complex '  &
            'Potentials')
END SUBROUTINE SO_Sector_Time_Propagator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
