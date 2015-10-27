!deck split_operator.f
!***begin prologue     split_operator
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           split-operator, propagation, RSP2, RSP4
!***author             schneider, b. i.(nsf)
!***source
!***purpose            time propagation using a split operator
!***                   based on Lie-Trotter-Suziki formulation.
!***                   A second and fourth order method have
!***                   been programmed.  the splitting is based
!***                   on the finite element DVR approach.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       split_operator
  Subroutine split_operator
  USE arnoldi_global_it
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  USE finite_element_propagator
  USE psi_h_psi
  USE h_on_vector
  USE pack_and_read_h
  IMPLICIT NONE
  REAL*4                                   :: secnds
  REAL*8                                   :: e_cur
  INTEGER                                  :: i, j
  INTEGER, DIMENSION(2)                    :: words
!
!                     Propagation
!
  CALL setup
!
!           Allocate the storage for the diagonal arrays
!
  ALLOCATE(v_tot(n3d),exp_diag(n3d))

  IF(typke == 'dvr'.OR.typke == 'packed') THEN
     IF (e_method == 'hamiltonian') THEN
!
!           Calculate Energy in propagation using Hamiltonian
!
         ALLOCATE(buf(spdim))
         DO i=1,spdim
            ALLOCATE( buf(i)%d(nphy(i)),              &
                      buf(i)%hbuf(nphy(i)*nphy(i)),   &
                      buf(i)%hibuf(2,nphy(i)*nphy(i)))
            CALL pack_h(i)
         END DO
     END IF
  END IF
!
!           Allocate needed storage for the one, two or three dimensional routines.
!
  IF(spdim == 1) THEN
     ALLOCATE(psi_1d(nphy(1)),v_scr_1d(nphy(1)))
     CALL so_prop(psi_1d,v_scr_1d)
  ELSE IF(spdim == 2 ) THEN
     ALLOCATE(psi_2d(nphy(2),nphy(1)),v_scr_2d(nphy(2),nphy(1)))
     maxdim=max(nphy(2),nphy(1))
     ALLOCATE( f_1(nphy(1)), f_2(maxdim*2,1) )
     CALL so_prop(psi_2d,v_scr_2d)
  ELSE IF(spdim == 3 ) THEN
     ALLOCATE(psi_3d(nphy(3),nphy(2),nphy(1)),                    &
              v_scr_3d(nphy(3),nphy(2),nphy(1)))
     maxdim=max(nphy(3),nphy(2),nphy(1))
     words(1)=max(nphy(2)*nphy(1),maxdim*2)
     words(2)=maxdim*maxdim*2
     ALLOCATE( f_1(nphy(1)), f_2(words(1),1), f_3(words(2),1,1) )
     CALL so_prop(psi_3d,v_scr_3d)
  ELSE
     CALL lnkerr('error in dimension')
  END IF
!
!           Deallocate all arrays that have been allocated for the entire
!           calculation.
! 
  DEALLOCATE(exp_diag,tim_pts)
  IF (e_method == 'exponential') THEN
     IF(typke == 'dvr'.OR.typke == 'packed') THEN
         ALLOCATE(buf(spdim))
         DO i=1,spdim
            ALLOCATE( buf(i)%d(nphy(i)), &
                      buf(i)%hbuf(nphy(i)*nphy(i)), &
                      buf(i)%hibuf(2,nphy(i)*nphy(i)))
            CALL pack_h(i)
         END DO
     END IF
     IF(spdim == 1) THEN
        call h_v_d(psi_1d,v_scr_1d,1)
        call check_energy(psi_1d,v_scr_1d,e_cur)
        DEALLOCATE(psi_1d,v_scr_1d)
     ELSE IF(spdim == 2 ) THEN
        call h_v_d(psi_2d,v_scr_2d,1)
        call check_energy(psi_2d,v_scr_2d,f_1,e_cur)
        DEALLOCATE(psi_2d,v_scr_2d)
        DEALLOCATE( f_1,f_2 )
     ELSE IF(spdim == 3 ) THEN
        call h_v_d(psi_3d,v_scr_3d,1)
        call check_energy(psi_3d,v_scr_3d,f_1,f_2,e_cur)
        DEALLOCATE(psi_3d,v_scr_3d)
        DEALLOCATE( f_1,f_2,f_3 )
     END IF
  END IF
  IF(typke /= 'fd') THEN
     DO i=1,spdim
        DEALLOCATE(grid(i)%pt,grid(i)%wt,grid(i)%f,grid(i)%ke,grid(i)%v)
     END DO
  ELSE
     DO i=1,spdim
        DEALLOCATE(grid(i)%pt,grid(i)%wt,grid(i)%ke,grid(i)%v)
     END DO
  END IF
  IF(diag) THEN
     DO i=1,spdim
        DEALLOCATE(grid(i)%eigv,grid(i)%eigvec,grid(i)%eigv_0,    &
                   grid(i)%eigvec_0)
     END DO
  END IF
  DO i=1,spdim
     DO j=1,num_reg(i)
        DEALLOCATE(mat_reg_d(j,i)%exp_t_mat)
     END DO
  END DO
  DEALLOCATE( grid,mat_reg_d,num_reg,nfun_reg )
  IF(typke == 'dvr'.OR.typke == 'packed') THEN
     DO i=1,spdim
        DEALLOCATE( buf(i)%d,                         &
                    buf(i)%hbuf,                      &
                    buf(i)%hibuf)
     END DO
     DEALLOCATE(buf)
  END IF
END SUBROUTINE split_operator

