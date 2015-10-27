!deck regional_matrices.f
!***begin prologue     regional_matrices
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           space, propagation, dvr, finite difference
!***author             schneider, b. i.(nsf)
!***source
!***purpose            the global FEDVR matrices are split into their
!***                   regional parts for ease of use in the propagation.
!***                   the regional kinetic energy matrix may be modified
!***                   to include any diagonal contribution from the one
!***                   body potential or any other diagnal part.  once this
!***                   done, the kinetic energy is diagonalized.  the 
!***                   eigenvectors and eigenvalues are used to form the
!***                   element propagators later in the code.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       space
  Subroutine Regional_matrices
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  USE regional_diagonalization 
  IMPLICIT NONE
  INTEGER                                  :: i, j, start, row_dim
  CHARACTER(LEN=8)                         :: mat_typ
!
!
  ALLOCATE(mat_reg(maxreg,spdim))
  IF(algorithm == 'split_operator') THEN
     DO  i=1,spdim
         IF(typke == 'dvr'.OR.typke == 'packed') THEN
            write(iout,1) i
         ELSE
            write(iout,2) i
         END IF
         WRITE(iout,3)
         maxmem_reg=0
         DO j=1,num_reg(i)
            write(iout,4) j, nfun_reg(j,i)
            ALLOCATE(mat_reg(j,i)%pt_d                             &
                    (nfun_reg(j,i)),                               &
                     mat_reg(j,i)%ke_mat_d                         &
                    (nfun_reg(j,i),nfun_reg(j,i)),                 &
                     mat_reg(j,i)%eigval_mat_d                     &
                    (nfun_reg(j,i)),                               &
                     mat_reg(j,i)%eigvec_mat_d                     &
                    (nfun_reg(j,i),nfun_reg(j,i)))
            maxmem_reg=max(maxmem_reg,nfun_reg(j,i))
         END DO
         ALLOCATE(scr_d(5*maxmem_reg))
         IF(typke == 'dvr'.OR.typke == 'packed') THEN
            mat_typ='full'
            row_dim=nphy(i)
         ELSE
            mat_typ='banded'
            row_dim=row(i)
         END IF
         IF(diag_mod /= 'none') THEN
            call modify_diag(grid(i)%ke,grid(i)%v,row_dim,nphy(i), &
                             diag_mod,mat_typ)
         END IF
         start = 1
         DO j=1,num_reg(i)
            call pt_reg(grid(i)%pt(start),mat_reg(j,i)%pt_d,       &
                        nfun_reg(j,i),nphy(i),j)
            start = start + nfun_reg(j,i) - 1
         END DO
!
!        here we actually form the FEDVR kinetic energy regional
!        matrices.
!
         start = 1
         IF(mat_typ=='full') THEN
            DO j=1,num_reg(i)
               call ke_reg_dvr(grid(i)%ke(start,start),            &
                               mat_reg(j,i)%ke_mat_d,              &
                               nfun_reg(j,i),nphy(i),j)
               start = start + nfun_reg(j,i) - 1
            END DO
         ELSE IF(mat_typ=='banded') THEN
            DO j=1,num_reg(i)
               call ke_reg_fd(grid(i)%ke(1,start),                 &
                              mat_reg(j,i)%ke_mat_d,               &
                              nfun_reg(j,i),row_dim,j)
               start = start + nfun_reg(j,i) - 1
            END DO
         ELSE
            CALL lnkerr('error')
         END IF  

         DO j=1,num_reg(i)
            write(iplot(1),*) 'real sector = ', j, 'size = ',npt(j)
            write(iplot(1),*) 'kinetic energy minus diagonals'
            write(iplot(1),*) mat_reg(j,i)%ke_mat_d
!
!           Diagonalize
!
            call regional_diag(mat_reg(j,i)%ke_mat_d,                   &
                               mat_reg(j,i)%eigval_mat_d,               &
                               mat_reg(j,i)%eigvec_mat_d,               &
                               nfun_reg(j,i),j)
         END DO
         DEALLOCATE(scr_d)
     END DO
  ELSE
     DO  i=1,spdim
         IF(typke == 'dvr'.OR.typke == 'packed') THEN
            write(iout,1) i
         ELSE
            write(iout,2) i
         END IF
         WRITE(iout,3)
         DO j=1,num_reg(i)
            write(iout,4) j, nfun_reg(j,i)
            ALLOCATE(mat_reg(j,i)%pt_d                             &
                    (nfun_reg(j,i)),                               &
                     mat_reg(j,i)%ke_mat_d                         &
                    (nfun_reg(j,i),nfun_reg(j,i)))
         END DO
         IF(typke == 'dvr'.OR.typke == 'packed') THEN
             mat_typ='full'
             row_dim=nphy(i)
         ELSE
             mat_typ='banded'
             row_dim=row(i)
         END IF
         IF(diag_mod /= 'none') THEN
            call modify_diag(grid(i)%ke,grid(i)%v,row_dim,nphy(i), &
                             diag_mod,mat_typ)
         END IF
         start = 1
         DO j=1,num_reg(i)
            call pt_reg(grid(i)%pt(start),mat_reg(j,i)%pt_d,       &
                        nfun_reg(j,i),nphy(i),j)
            start = start + nfun_reg(j,i) - 1
         END DO
         start = 1
         IF(mat_typ=='full') THEN
            DO j=1,num_reg(i)
               call ke_reg_dvr(grid(i)%ke(start,start),            &
                               mat_reg(j,i)%ke_mat_d,              &
                               nfun_reg(j,i),nphy(i),j)
               start = start + nfun_reg(j,i) - 1
            END DO
         ELSE IF(mat_typ=='banded') THEN
            DO j=1,num_reg(i)
               call ke_reg_fd(grid(i)%ke(1,start),                 &
                              mat_reg(j,i)%ke_mat_d,               &
                              nfun_reg(j,i),row_dim,j)
               start = start + nfun_reg(j,i) - 1
            END DO
         ELSE
            CALL lnkerr('error')
         END IF  
!
         DO j=1,num_reg(i) 
            write(iplot(1),*) 'sector = ', j, 'size = ',npt(j)
            write(iplot(1),*) 'kinetic energy minus diagonals'
            write(iplot(1),*) mat_reg(j,i)%ke_mat_d
         END DO
     END DO
  END IF
1 FORMAT(/15x,'Regional Information for DVR Basis,'                &
              ' Variable = ',i2)
2 FORMAT(/15x,'Regional Information for FD Basis,'                 &
                 ' Variable = ',i2)
3 FORMAT(/,17x,'Region',5x,'Number of Functions')
4 FORMAT(16x,i4,14x,i5)
END Subroutine regional_matrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
