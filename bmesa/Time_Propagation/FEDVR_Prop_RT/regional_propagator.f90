!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      MODULE regional_propagator
!**begin prologue     regional_propagator
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       regional_propagator
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck reg_prop_d
!***begin prologue     reg_prop_d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            This routine performs the basic splitting operation on a
!***                   Hamiltonian containing multiple subregions.  We label the
!***                   regions (1,3,5,....) as odd and the regions (2,4,6,....) as even.
!***                   H_a is the matrix which is the sum of the odd regions and H_b
!***                   the sum of the even regions.
!
!***description        The time propagators for each of the subregions is constructed
!***                   from the eigenvalues and eigenvectors of the regional kinetic
!***                   energy matrices.
!***references
!***routines called
!***end prologue       reg_prop_d
!
  SUBROUTINE reg_prop_d(first_region,last_region,q,pointer)
  USE dvr_global
  USE dvrprop_global_rt
  USE regional_umat
  IMPLICIT NONE
  INTEGER                                  :: first_region, last_region
  INTEGER                                  :: q, pointer
  DO i=first_region,last_region,2
          call umat_reg(mat_reg_d(i,q)%eigval_mat_d,                    &
                        mat_reg_d(i,q)%eigvec_mat_d,                    &
                        mat_reg_d(i,q)%eigvec_mat_d,                    &
                        mat_reg_d(i,q)%cosine_t_mat(:,:,pointer),       &
                        mat_reg_d(i,q)%sine_t_mat(:,:,pointer),         &
                        si_d,ci_d,tau,nfun_reg(i,q))
  END DO
END SUBROUTINE reg_prop_d
!deck reg_prop_z
!***begin prologue     reg_prop_z
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            This routine performs the basic splitting operation on a
!***                   Hamiltonian containing multiple subregions.  We label the
!***                   regions (1,3,5,....) as odd and the regions (2,4,6,....) as even.
!***                   H_a is the matrix which is the sum of the odd regions and H_b
!***                   the sum of the even regions.
!
!***description        The time propagators for each of the subregions is constructed
!***                   from the eigenvalues and eigenvectors of the regional kinetic
!***                   energy matrices.
!***references
!***routines called
!***end prologue       reg_prop_z
!
  SUBROUTINE reg_prop_z(first_region,last_region,q,pointer)
  USE dvr_global
  USE dvrprop_global_rt
  USE regional_umat
  IMPLICIT NONE
  INTEGER                                  :: first_region, last_region
  INTEGER                                  :: q, pointer
  DO i=first_region,last_region,2
          call umat_reg(mat_reg_z(i,q)%eigval_mat_z,                &
                        mat_reg_z(i,q)%eigvec_mat_z_r,              &
                        mat_reg_z(i,q)%eigvec_mat_z_l,              &
                        mat_reg_z(i,q)%cosine_t_mat(:,:,pointer),   &
                        mat_reg_z(i,q)%sine_t_mat(:,:,pointer),     &
                        si_z,ci_z,tau,nfun_reg(i,q))
  END DO
END SUBROUTINE reg_prop_z
END MODULE regional_propagator















