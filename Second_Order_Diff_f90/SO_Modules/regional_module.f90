!***********************************************************************
! regional_module
!**begin prologue     regional_module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Contains all of the subroutines to form,
!***                  diagonalize and construct the regional matrices.
!***                  Explicit interfaces are used to allow
!***                  a transparent use of generic subroutines which work
!***                  for both real and complex vectors.  This feature
!***                  permits a single code to be used for both real and
!***                  imaginary time propagation.
!***description       See subroutines
!**references
!**modules needed     See USE statements below
!**end prologue       regional_module
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      MODULE regional_module
                      USE dvr_global
                      USE dvr_shared
                      USE dvrprop_global
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!***                   is done, the kinetic energy is diagonalized.  the 
!***                   eigenvectors and eigenvalues are used to form the
!***                   element propagators later in the code.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       space
  Subroutine Regional_matrices
  IMPLICIT NONE
  INTEGER                                  :: i, j, start, row_dim
  CHARACTER(LEN=8)                         :: mat_typ
!
!
  mat_typ='full'
  ALLOCATE(mat_reg(maxreg,spdim))
  DO  i=1,spdim
      write(iout,1) i
      WRITE(iout,2)
      DO j=1,n_reg_real(i)
         write(iout,3) j, nfun_reg(j,i), 'real'
         ALLOCATE(mat_reg(j,i)%pt_d                             &
                 (nfun_reg(j,i)),                               &
                  mat_reg(j,i)%ke_mat_d                         &
                 (nfun_reg(j,i),nfun_reg(j,i)))
      END DO
      row_dim=nphy(i)
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
      DO j=1,num_reg(i)
           call ke_reg_dvr(grid(i)%ke(start,start),              &
                           mat_reg(j,i)%ke_mat_d,                &
                           nfun_reg(j,i),nphy(i),j)
           start = start + nfun_reg(j,i) - 1
      END DO
  END  DO
1 FORMAT(/15x,'Regional Information for DVR Basis,'                &
              ' Variable = ',i2)
2 FORMAT(/,17x,'Region',5x,'Number of Functions',5x,'Type')
3 FORMAT(16x,i4,14x,i5,9x,a8)
END Subroutine regional_matrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE regional_module
