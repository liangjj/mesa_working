!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              MODULE exp_off_diagonal
!deck exp_off_diagonal
!***begin prologue     exp_off_diagonal     
!***date written       040607   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            compute the effect of the split operator
!                      off-diagonal, exponential propagator on a vector.
!
!***description        The operation consists of applying,
!
!     exp(-H_odd*t/(2*hbar)) * exp(-H_even*t/(hbar)) * exp(-H_odd*t/(2*hbar))
!
!                      to the vector.  H_even and H_odd are separately block
!                      diagonal but the blocks overlap.  Thus the operations must
!                      alternate their input vectors.
!
!                     H = H(1,1) + H(2,2) + H(3,3)
!
!     Note that vectors are stored as V(nphy(3),nphy(2),nphy(1)).  Thus matrix
!     vector multiply in two and three dimensions must be done carefully for
!     nphy(2) and nphy(3).  We have( repeated indices summed over), in 3D,
!    
! V(3,2,1) = [ H(1,1) * V(3,2,1) + H(2,2) * V(3,2,1) + H(3,3) * V(3,2,1) ]  
!                                   =
!            [ V(3,2,1) * H(1,1) + V(3,2,1) * H(2,2) + H(3,3) * V(3,2,1) ]
!
! Thus the first mutiply may be done as a matrix V(3*2,1) on H(1,1), 
! the second as V(3,2,1) * H(2,2) with an outer loop over index 1 and 
! the third as a simple matrix vector multiply, H(3,3) * V(3,2*1).
!
!and in 2D,
!
! V(2,1)   = [ H(1,1) * V(2,1) + H(2,2) * V(2,1)  ]  
!                                   =
!            [ V(2,1) * H(1,1) + H(2,2) * V(2,1) ]
!
!
!
!***references
!***routines called    dvr_exp_off_diag_m_v, fd_exp_off_diag_m_v
!***                   dvr_exp_off_diag_v_m_2_d, fd_exp_off_diag_v_m_2_d
!***                   dvr_exp_off_diag_v_m_3_d, fd_exp_off_diag_v_m_3_d
!
!***end prologue       exp_off_diagonal
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              INTERFACE exp_off_diag
                    MODULE PROCEDURE exp_off_diag_1d_d,       &
                                     exp_off_diag_2d_d,       &
                                     exp_off_diag_3d_d
                          END INTERFACE exp_off_diag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              CONTAINS 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE exp_off_diag_1d_d(wave_function,scratch_vector)
     USE dvrprop_global_it
     USE dvr_shared
     USE dvr_global
     IMPLICIT NONE
     REAL*8, DIMENSION(nphy(1))                   :: wave_function
     REAL*8, DIMENSION(nphy(1))                   :: scratch_vector
     call exp_off_diag_m_v(wave_function,          &
                           scratch_vector,         &
                           nphy(1),1,1)
  END SUBROUTINE exp_off_diag_1d_d
  SUBROUTINE exp_off_diag_2d_d(wave_function,scratch_vector)
     USE dvrprop_global_it
     USE dvr_shared
     USE dvr_global
     IMPLICIT NONE
     REAL*8, DIMENSION(nphy(2),nphy(1))           :: wave_function
     REAL*8, DIMENSION(nphy(2),nphy(1))           :: scratch_vector
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!            This code is for the general problem involving nphy(2),
!            that is H(2,2) * V(2,1) 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      call exp_off_diag_m_v(wave_function,          &
                            scratch_vector,         &
                            nphy(2),nphy(1),2)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!            This code is for the general problem involving nphy(1),
!            V(2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      call exp_off_diag_v_m_2_d(wave_function,      &
                                scratch_vector,     &
                                nphy(2),nphy(1),1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END SUBROUTINE exp_off_diag_2d_d
  SUBROUTINE exp_off_diag_3d_d(wave_function,scratch_vector)
     USE dvrprop_global_it
     USE dvr_shared
     USE dvr_global
     USE plot_wavefunction
     IMPLICIT NONE
     REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))   :: wave_function
     REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))   :: scratch_vector

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            This code is for the general problem involving nphy(3),
!            H(3,3) * V(3,2,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      call exp_off_diag_m_v(wave_function,         &
                            scratch_vector,        &
                            nphy(3),nphy(2)*nphy(1),3)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general case involving the matrix multiply 
!        involving nphy(2), that is, V(3,2,1) * H(2,2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call exp_off_diag_v_m_3_d                    &
                        (wave_function,            &
                         scratch_vector,           &
                         nphy(3),nphy(2),nphy(1),2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general problem involving the matrix multiply 
!        for nphy(1), that is, V(3,2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call exp_off_diag_v_m_2_d                   &
                       (wave_function,            &
                        scratch_vector,           &
                        nphy(3)*nphy(2),nphy(1),1)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END SUBROUTINE exp_off_diag_3d_d
!
!deck exp_off_diag_m_v
!***begin prologue     exp_off_diag_m_v     
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
!                 V = exp_t_mat * V
!
!------------------------------------------------------------------------------------
!***                   where exp_t_mat is the regional
!***                   propagator for a FEDVR or FD Hamiltonian.
!
!***                   The parameter nj is a dummy and can take on values
!***                   consistent with any dimensional problem.  In a 1D
!***                   case nj=1, in 2D nj=nx and in 3D_D nj=ny*nx.

!***references
!***routines called    v_m_v_gen, v_m_v_2, v_m_v_3
!
!***end prologue       exp_off_diag_m_v                     
!
  SUBROUTINE exp_off_diag_m_v(v,           &
                              v_scr,       &
                              ni,nj,index)
  USE io
  USE dvr_global,                     ONLY    : num_reg, nfun_reg
  USE dvrprop_global_it,              ONLY    : mat_reg_d, prop_point
  USE v_m_v
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, index
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  INTEGER                                  :: i, ntrips, trips
  INTEGER, DIMENSION(3)                    :: locate, begin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
!  The sequence is called three times.  First for the odd matrices, then
!  for the even matrices and then again for the odd matrices.
! 
  locate(1) = 1
  locate(2) = nfun_reg(1,index)
  locate(3) = 1
  begin(1) = 1
  begin(2) = 2
  begin(3) = 1
  ntrips=3
  IF(num_reg(index) == 1 ) THEN
     ntrips=1
  END IF
  DO trips=1,ntrips
     DO i = begin(trips), num_reg(index), 2
!
        IF ( nfun_reg(i,index) > 7 ) then
!
!                       General Code
!
           CALL v_m_v_gen(v(locate(trips):ni,1:nj),                           &
                          v_scr(locate(trips):ni,1:nj),                       &
                          mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),       &
                          ni,nj,nfun_reg(i,index))
        ELSE IF (nfun_reg(i,index) == 2) then
! 
!                       Special case for Two by Two
!
           CALL v_m_v_2(v(locate(trips):ni,1:nj),                             &
                        v_scr(locate(trips):ni,1:nj),                         &
                        mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),         &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Special case for Three by Three
!
           CALL v_m_v_3(v(locate(trips):ni,1:nj),                             &
                        v_scr(locate(trips):ni,1:nj),                         &
                        mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),         &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Special case for Four by Four
!
           CALL v_m_v_4(v(locate(trips):ni,1:nj),                             &
                        v_scr(locate(trips):ni,1:nj),                         &
                        mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),         &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Special case for Five by Five
!
           CALL v_m_v_5(v(locate(trips):ni,1:nj),                             &
                        v_scr(locate(trips):ni,1:nj),                         &
                        mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),         &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Special case for Six by Six
!
           CALL v_m_v_6(v(locate(trips):ni,1:nj),                             &
                        v_scr(locate(trips):ni,1:nj),                         &
                        mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),         &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Special case for Seven by Seven
!
           CALL v_m_v_7(v(locate(trips):ni,1:nj),                             &
                        v_scr(locate(trips):ni,1:nj),                         &
                        mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),         &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Special case for Eight by Eight
!
           CALL v_m_v_8(v(locate(trips):ni,1:nj),                             &
                        v_scr(locate(trips):ni,1:nj),                         &
                        mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),         &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Special case for Nine by Nine
!
           CALL v_m_v_9(v(locate(trips):ni,1:nj),                             &
                        v_scr(locate(trips):ni,1:nj),                         &
                        mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),         &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Special case for Ten by Ten
!
           CALL v_m_v_10(v(locate(trips):ni,1:nj),                            &
                         v_scr(locate(trips):ni,1:nj),                        &
                         mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),        &
                         ni,nj)
        END IF
        IF (i /= num_reg(index) ) then
            locate(trips) = locate(trips) + nfun_reg(i,index)                 &
                                          +                                   &
                            nfun_reg(i+1,index)                               &
                                          -                                   &
                                          2
        END IF
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE exp_off_diag_m_v
!
!

!*deck exp_off_diag_v_m_2_d
!***begin prologue     exp_off_diag_v_m_2_d     
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
!                 V = V * exp_t_mat 
!
!------------------------------------------------------------------------------------
!***                   where exp_t_mat is the regional
!***                   propagator for a FEDVR or FD Hamiltonian. 
!***                   parts.  
!***                   The routine is not needed in 1D problems. 
!***                   In 2D, the ni index runs over the number
!                      of points in the y coordinate while in 3D_D,
!***                   it runs over the product nz*ny of the number 
!***                   of points in the z and y coordinates.
!
!***references
!***routines called    v_v_m_gen, v_v_m_2, v_v_m_3
!***end prologue       exp_off_diag_v_m_2_d
!
  SUBROUTINE exp_off_diag_v_m_2_d(v,       &
                                  v_scr,   &
                                  ni,nj,index)
  USE io
  USE dvr_global,                ONLY    : num_reg, nfun_reg
  USE dvrprop_global_it,         ONLY    : mat_reg_d, prop_point
  USE v_v_m
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, index, ntrips
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  INTEGER                                  :: i, trips
  INTEGER, DIMENSION(3)                    :: locate, begin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  locate(1) = 1
  locate(2) = nfun_reg(1,index)
  locate(3) = 1
  begin(1) = 1
  begin(2) = 2
  begin(3) = 1
  ntrips = 3
  IF(num_reg(index) == 1 ) THEN
     ntrips=1
  END IF
!
  DO trips = 1, ntrips
     DO i=begin(trips), num_reg(index), 2
        IF ( nfun_reg(i,index) > 7 ) then
!
!                        General Case
!
           CALL v_v_m_gen(v(1:ni,locate(trips):nj),                           &
                          v_scr(1:ni,locate(trips):nj),                       &
                          mat_reg_d(i,index)%exp_t_mat(:,:,prop_point), &
                          ni,nj,nfun_reg(i,index)) 
!
        ELSE IF (nfun_reg(i,index) == 2) then
!
!                       Two by Two
! 
           CALL v_v_m_2(v(1:ni,locate(trips):nj),                             &
                        v_scr(1:ni,locate(trips):nj),                         &
                        mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),   &
                        ni,nj)
       ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Three by Three
!
           CALL v_v_m_3(v(1:ni,locate(trips):nj),                             &
                        v_scr(1:ni,locate(trips):nj),                         &
                        mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),   &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Four by Four
!
           CALL v_v_m_4(v(1:ni,locate(trips):nj),                             &
                        v_scr(1:ni,locate(trips):nj),                         &
                        mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),   &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Five by Five
!
           CALL v_v_m_5(v(1:ni,locate(trips):nj),                             &
                        v_scr(1:ni,locate(trips):nj),                         &
                        mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),   &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Six by Six
!
           CALL v_v_m_6(v(1:ni,locate(trips):nj),                             &
                        v_scr(1:ni,locate(trips):nj),                         &
                        mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),   &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Seven by Seven
!
           CALL v_v_m_7(v(1:ni,locate(trips):nj),                             &
                        v_scr(1:ni,locate(trips):nj),                         &
                        mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),   &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Eight by Eight
!
           CALL v_v_m_8(v(1:ni,locate(trips):nj),                             &
                        v_scr(1:ni,locate(trips):nj),                         &
                        mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),   &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Nine by Nine
!
           CALL v_v_m_9(v(1:ni,locate(trips):nj),                             &
                        v_scr(1:ni,locate(trips):nj),                         &
                        mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),   &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Ten by Ten
!
           CALL v_v_m_10(v(1:ni,locate(trips):nj),                            &
                         v_scr(1:ni,locate(trips):nj),                        &
                         mat_reg_d(i,index)%exp_t_mat(:,:,prop_point),  &
                         ni,nj)
        END IF
        IF (i /= num_reg(index) ) then
            locate(trips) = locate(trips) + nfun_reg(i,index)           &
                                          +                             &
                            nfun_reg(i+1,index)                         &
                                          -                             &
                                          2
        END IF
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE exp_off_diag_v_m_2_d
!
!
!*deck exp_off_diag_v_m_3_d
!***begin prologue     exp_off_diag_v_m_3_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine simply calls
!***                   exp_off_diag_v_m_2_d in a loop over a dummy
!***                   index, k
!
!***references
!***routines called
!***end prologue       exp_off_diag_v_m_3_d
!
  SUBROUTINE exp_off_diag_v_m_3_d(v,      &
                                  v_scr,  &
                                 ni,nj,nk,index)
  IMPLICIT NONE
  INTEGER                             :: ni, nj, nk, index
  REAL*8, DIMENSION(ni,nj,nk)         :: v
  REAL*8, DIMENSION(ni,nj,nk)         :: v_scr
  INTEGER                             :: k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  DO k=1,nk
     CALL exp_off_diag_v_m_2_d(v(1,1,k),          &
                               v_scr(1,1,k),      &
                               ni,nj,index) 
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE exp_off_diag_v_m_3_d
END MODULE exp_off_diagonal

