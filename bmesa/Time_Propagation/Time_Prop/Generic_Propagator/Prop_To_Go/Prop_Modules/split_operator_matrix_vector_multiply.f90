!***********************************************************************
!split_operator_matrix_vector_multiply_module
!**begin prologue     split_operator_matrix_vector_multiply_module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Contains all of the major subroutines to perform
!***                  the matrix vector operation in the SO
!***                  algorithm.  Explicit interfaces are used to allow
!***                  a transparent use of generic subroutines which work
!***                  for both real and complex vectors.  This feature
!***                  permits a single code to be used for both real and
!***                  imaginary time propagation.
!***description       Given a starting vector, the wavefunction is propagated 
!***                  until a certain convergence criterion is reached.
!***                  In imaginary time, the wavefunction is propagated 
!***                  until the eigenvalue has converged to a given accuracy. 
!***                  For real time, one just proceeds from step to step. 
!**references
!**modules needed     See USE statements below
!**end prologue       split_operator_matrix_vector_multiply_module
!***********************************************************************
                     MODULE split_operator_matrix_vector_multiply_module
                     USE arnoldi_global
                     USE dvrprop_global
                     USE dvr_global
                     USE dvr_shared
                     USE regional_module
                     USE normalize_module
                     USE spatial_wavefunction_module
                     USE initial_state_module
                     USE moment_module
                     USE auto_correlation_module
                     USE non_linear_potential_module
                     USE matrix_vector_multiply_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     INTERFACE so_exp_diagonal_multiply
            MODULE PROCEDURE so_exp_diagonal_mul_d,                     &
                             so_exp_diagonal_mul_z
                 END INTERFACE so_exp_diagonal_multiply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     INTERFACE so_exp_off_diagonal_multiply
            MODULE PROCEDURE so_exp_off_diagonal_mul_d,         &
                             so_exp_off_diagonal_mul_z      
                 END INTERFACE so_exp_off_diagonal_multiply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck so_exp_diagonal_mul_d
!***begin prologue     so_exp_diagonal_mul_d
!***date written       040707   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            diagonal scaling of propagated vector by exponential.
!***
!***description        
!                      
!                      
!***references
!***routines called    
!***end prologue       so_exp_diagonal_mul_d
!
  SUBROUTINE so_exp_diagonal_mul_d(vector_in_out)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                         :: vector_in_out
  vector_in_out(:) = exp_diag_d(:) * vector_in_out(:)  
END SUBROUTINE so_exp_diagonal_mul_d
!***********************************************************************
!***********************************************************************
!deck so_exp_diagonal_mul_z
!***begin prologue     exp_diagonal_mul_z
!***date written       040707   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            diagonal scaling of propagated vector
!***
!***description        
!                      
!                      
!***references
!***routines called    
!***end prologue       so_exp_diagonal_mul_z
!
  SUBROUTINE so_exp_diagonal_mul_z(vector_in_out)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)                       :: vector_in_out
  vector_in_out(:) = exp_diag_z(:) * vector_in_out(:)  
END SUBROUTINE so_exp_diagonal_mul_z
!***********************************************************************
!***********************************************************************
!deck so_exp_off_diagonal_mul
!***begin prologue     so_exp_off_diagonal_mul     
!***date written       040607   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            compute the effect of the split operator
!                      off-diagonal, exponential propagator on a vector.
!
!***description        These routine perform the complex operation of
!***                   applying the exponential of a FEDVR Hamiltonian
!***                   on a vector.  The operation consists of applying,
!
!     exp(-iH_odd*t/(2*hbar)) * exp(-iH_even*t/(hbar)) * exp(-iH_odd*t/(2*hbar))
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
! In 2D, we have,
!
! V(2,1)   = [ H(1,1) * V(2,1) + H(2,2) * V(2,1)  ]  
!                                   =
!            [ V(2,1) * H(1,1) + H(2,2) * V(2,1) ]
!
!
!
!***references
!***routines called    
!***                  
!
!***end prologue       so_exp_off_diagonal_mul
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE so_exp_off_diagonal_mul_d(wave_function,scratch_vector)
     IMPLICIT NONE
     REAL*8, DIMENSION(:)                   :: wave_function
     REAL*8, DIMENSION(:)                   :: scratch_vector
     IF (spdim == 1) THEN
         CALL exp_off_diag_m_v_d(wave_function,                      &
                                 scratch_vector,                     &
                                 nphy(1),1,1)
     ELSE IF ( spdim == 2) THEN
         CALL so_exp_off_diagonal_mul_2d_d(wave_function,            &
                                           scratch_vector)
     ELSE IF( spdim == 3) THEN
         CALL so_exp_off_diagonal_mul_3d_d(wave_function,            &
                                           scratch_vector)
     END IF
  END SUBROUTINE so_exp_off_diagonal_mul_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE so_exp_off_diagonal_mul_z(wave_function,scratch_vector)
     IMPLICIT NONE
     COMPLEX*16, DIMENSION(:)                   :: wave_function
     COMPLEX*16, DIMENSION(:)                   :: scratch_vector
     IF (spdim == 1) THEN
         CALL exp_off_diag_m_v_z(wave_function,                      &
                                 scratch_vector,                     &
                                 nphy(1),1,1)
     ELSE IF ( spdim == 2) THEN
         CALL so_exp_off_diagonal_mul_2d_z(wave_function,            &
                                           scratch_vector)
     ELSE IF( spdim == 3) THEN
         CALL so_exp_off_diagonal_mul_3d_z(wave_function,            &
                                           scratch_vector)
     END IF
  END SUBROUTINE so_exp_off_diagonal_mul_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE so_exp_off_diagonal_mul_2d_d(wave_function,scratch_vector)
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
      call exp_off_diag_m_v_d(wave_function,                     &
                              scratch_vector,                    &
                              nphy(2),nphy(1),2)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!            This code is for the general problem involving nphy(1),
!            V(2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      call exp_off_diag_v_m_2_d(wave_function,                   &
                                scratch_vector,                  &
                                nphy(2),nphy(1),1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END SUBROUTINE so_exp_off_diagonal_mul_2d_d
!***********************************************************************
  SUBROUTINE so_exp_off_diagonal_mul_3d_d(wave_function,scratch_vector)
     USE dvrprop_global
     USE dvr_shared
     USE dvr_global
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
      call exp_off_diag_m_v_d(wave_function,                     &
                              scratch_vector,                    &
                              nphy(3),nphy(2)*nphy(1),3)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general case involving the matrix multiply 
!        involving nphy(2), that is, V(3,2,1) * H(2,2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call exp_off_diag_v_m_3_d                                  &
                             (wave_function,                     &
                              scratch_vector,                    &
                              nphy(3),nphy(2),nphy(1),2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general problem involving the matrix multiply 
!        for nphy(1), that is, V(3,2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call exp_off_diag_v_m_2_d                                  &
                             (wave_function,                     &  
                              scratch_vector,                    &
                              nphy(3)*nphy(2),nphy(1),1)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END SUBROUTINE so_exp_off_diagonal_mul_3d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck exp_off_diag_m_v_d
!***begin prologue     exp_off_diag_m_v_d     
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
!***end prologue       exp_off_diag_m_v_d                     
!
  SUBROUTINE exp_off_diag_m_v_d(v,                                             &
                                v_scr,                                         &
                                ni,nj,index)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, index
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  INTEGER                                  :: i, ntrips, trips
  INTEGER                                  :: first, last
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
!  write(iout,*) ni, nj, index
  DO trips=1,ntrips
     first = locate(trips)
     DO i = begin(trips), num_reg(index), 2
!         write(iout,*) trips, begin(trips), num_reg(index)
!         write(iout,*) i, nfun_reg(i,index), locate(trips), prop_point
!
        IF ( nfun_reg(i,index) > 10 ) then
!
!                       General Code
!
           last = first + nfun_reg(i,index) - 1
           CALL v_so_mat_v_gen(v(first:last,1:nj),                               &
                               v_scr(first:last,1:nj),                            &
                               mat_reg(i,index)%exp_t_mat_d(:,:,prop_point))
        ELSE IF (nfun_reg(i,index) == 2) then
! 
!                       Special case for Two by Two
!
!           write(iout,*) 'Input Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
           last = first + 1
           CALL v_so_mat_v_2(v(first:last,1:nj),                                &
                             v_scr(first:last,1:nj),                            &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point))
!           write(iout,*) 'Output Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
        ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Special case for Three by Three
!
!           write(iout,*) 'Input Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
           last = first + 2
           CALL v_so_mat_v_3(v(first:last,1:nj),                               &
                             v_scr(first:last,1:nj),                           &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point))
!           write(iout,*) 'Output Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
        ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Special case for Four by Four
!
           last = first + 3              
           CALL v_so_mat_v_4(v(first:last,1:nj),                                &
                             v_scr(first:last,1:nj),                            &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point))
        ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Special case for Five by Five
!
           last = first + 4
           CALL v_so_mat_v_5(v(first:last,1:nj),                                &
                             v_scr(first:last,1:nj),                            &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point))
        ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Special case for Six by Six
!
!           write(iout,*) 'Input Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
           last = first + 5
           CALL v_so_mat_v_6(v(first:last,1:nj),                                &
                             v_scr(first:last,1:nj),                            &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point))
!           write(iout,*) 'Output Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
        ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Special case for Seven by Seven
!
           last = first + 6
           CALL v_so_mat_v_7(v(first:last,1:nj),                                &
                             v_scr(first:last,1:nj),                            &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point))
        ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Special case for Eight by Eight
!
           last = first + 7
           CALL v_so_mat_v_8(v(first:last,1:nj),                                &
                             v_scr(first:last,1:nj),                            &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point))
        ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Special case for Nine by Nine
!
           last = first + 8
           CALL v_so_mat_v_9(v(first:last,1:nj),                                &
                             v_scr(first:last,1:nj),                            &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point))
        ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Special case for Ten by Ten
!
           last = first + 9
           CALL v_so_mat_v_10(v(first:last,1:nj),                               &
                              v_scr(first:last,1:nj),                           &
                              mat_reg(i,index)%exp_t_mat_d(:,:,prop_point))
        END IF
        IF (i /= num_reg(index) ) then
            first = first + nfun_reg(i,index)                                   &
                          +                                                     &
                            nfun_reg(i+1,index)                                 &
                          -                                                     &
                            2
        END IF
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE exp_off_diag_m_v_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  SUBROUTINE exp_off_diag_v_m_2_d(v,                                         &
                                  v_scr,                                     &
                                  ni,nj,index)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, index, ntrips
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  INTEGER                                  :: i, trips
  INTEGER, DIMENSION(3)                    :: locate, begin
  INTEGER                                  :: first, last
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
     first = locate(trips)
     DO i=begin(trips), num_reg(index), 2
        IF ( nfun_reg(i,index) > 10 ) then
!
!                        General Case
!
           last = first + nfun_reg(i,index) - 1
           CALL v_v_so_mat_gen(v(1:ni,first:last),                             &
                               v_scr(1:ni,first:last),                         &
                               mat_reg(i,index)%exp_t_mat_d(:,:,prop_point)) 
!
        ELSE IF (nfun_reg(i,index) == 2) then
!
!                       Two by Two
! 
           last= first + 1
           CALL v_v_so_mat_2(v(1:ni,first:last),                              &
                             v_scr(i:ni,first:last),                          &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point)) 
       ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Three by Three
!
           last = first + 2
           CALL v_v_so_mat_3(v(1:ni,first:last),                              &
                             v_scr(1:ni,first:last),                          &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point)) 
        ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Four by Four
!
           last= first + 3
           CALL v_v_so_mat_4(v(1:ni,first:last),                              &
                             v_scr(1:ni,first:last),                          &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point)) 
        ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Five by Five
!
           last= first + 4
           CALL v_v_so_mat_5(v(1:ni,first:last),                              &
                             v_scr(1:ni,first:last),                          &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point)) 
        ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Six by Six
!
           last= first + 5
           CALL v_v_so_mat_6(v(1:ni,first:last),                               &
                             v_scr(1:ni,first:last),                           &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point)) 
        ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Seven by Seven
!
           last= first + 6
           CALL v_v_so_mat_7(v(1:ni,first:last),                               &
                             v_scr(1:ni,first:last),                           &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point)) 
        ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Eight by Eight
!
           last= first + 7
           CALL v_v_so_mat_8(v(1:ni,first:last),                               &
                             v_scr(1:ni,first:last),                           &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point)) 
        ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Nine by Nine
!
           last= first + 8
           CALL v_v_so_mat_9(v(1:ni,first:last),                               &
                             v_scr(1:ni,first:last),                           &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point)) 
        ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Ten by Ten
!
           last= first + 9
           CALL v_v_so_mat_10(v(1:ni,first:last),                              &
                              v_scr(1:ni,first:last),                          &
                              mat_reg(i,index)%exp_t_mat_d(:,:,prop_point)) 
        END IF
        IF (i /= num_reg(index) ) then
            locate(trips) = locate(trips) + nfun_reg(i,index)                     &
                                          +                                       &
                            nfun_reg(i+1,index)                                   &
                                          -                                       &
                                          2
        END IF
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE exp_off_diag_v_m_2_d
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
!***                   exp_off_diag_v_m_2 in a loop over a dummy
!***                   index, k
!
!***references
!***routines called
!***end prologue       exp_off_diag_v_m_3_d
!
  SUBROUTINE exp_off_diag_v_m_3_d(v,                                         &
                                  v_scr,                                     &
                                  ni,nj,nk,index)
  IMPLICIT NONE
  INTEGER                             :: ni, nj, nk, index
  REAL*8, DIMENSION(ni,nj,nk)         :: v
  REAL*8, DIMENSION(ni,nj,nk)         :: v_scr
  INTEGER                             :: k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  DO k=1,nk
     CALL exp_off_diag_v_m_2_d(v(:,:,k),                                     &
                               v_scr(:,:,k),                                 &
                               ni,nj,index) 
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE exp_off_diag_v_m_3_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck so_exp_off_diagonal_mul_z
!***begin prologue     so_exp_off_diagonal_mul_z     
!***date written       040607   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            compute the effect of the split operator
!                      off-diagonal, exponential propagator on a vector.
!
!***description
!***references
!***routines called    
!***                  
!
!***end prologue       so_exp_off_diagonal_mul_z
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE so_exp_off_diagonal_mul_1d_z(wave_function,scratch_vector)
     IMPLICIT NONE
     COMPLEX*16, DIMENSION(nphy(1))                   :: wave_function
     COMPLEX*16, DIMENSION(nphy(1))                   :: scratch_vector
     call exp_off_diag_m_v_z(wave_function,                      &
                             scratch_vector,                     &
                             nphy(1),1,1)
  END SUBROUTINE so_exp_off_diagonal_mul_1d_z
  SUBROUTINE so_exp_off_diagonal_mul_2d_z(wave_function,scratch_vector)
     IMPLICIT NONE
     COMPLEX*16, DIMENSION(nphy(2),nphy(1))           :: wave_function
     COMPLEX*16, DIMENSION(nphy(2),nphy(1))           :: scratch_vector
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!            This code is for the general problem involving nphy(2),
!            that is H(2,2) * V(2,1) 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      call exp_off_diag_m_v_z(wave_function,                     &
                              scratch_vector,                    &
                              nphy(2),nphy(1),2)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!            This code is for the general problem involving nphy(1),
!            V(2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      call exp_off_diag_v_m_2_z(wave_function,                   &
                                scratch_vector,                  &
                                nphy(2),nphy(1),1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END SUBROUTINE so_exp_off_diagonal_mul_2d_z
  SUBROUTINE so_exp_off_diagonal_mul_3d_z(wave_function,scratch_vector)
     IMPLICIT NONE
     COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1))   :: wave_function
     COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1))   :: scratch_vector

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            This code is for the general problem involving nphy(3),
!            H(3,3) * V(3,2,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      call exp_off_diag_m_v_z(wave_function,                     &
                              scratch_vector,                    &
                              nphy(3),nphy(2)*nphy(1),3)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general case involving the matrix multiply 
!        involving nphy(2), that is, V(3,2,1) * H(2,2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call exp_off_diag_v_m_3_z                                  &
                             (wave_function,                     &
                              scratch_vector,                    &
                              nphy(3),nphy(2),nphy(1),2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general problem involving the matrix multiply 
!        for nphy(1), that is, V(3,2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call exp_off_diag_v_m_2_z                                  &
                             (wave_function,                     &
                              scratch_vector,                    &
                              nphy(3)*nphy(2),nphy(1),1)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END SUBROUTINE so_exp_off_diagonal_mul_3d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck exp_off_diag_m_v_z
!***begin prologue     exp_off_diag_m_v_z     
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
!***end prologue       exp_off_diag_m_v_z                     
!
  SUBROUTINE exp_off_diag_m_v_z(v,                                             &
                                v_scr,                                         &
                                ni,nj,index)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, index
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  INTEGER                                  :: i, ntrips, trips
  INTEGER                                  :: first, last
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
  write(iout,*) ni, nj, index
  DO trips=1,ntrips
     first = locate(trips)
     DO i = begin(trips), num_reg(index), 2
!         write(iout,*) trips, begin(trips), num_reg(index)
!         write(iout,*) i, nfun_reg(i,index), locate(trips), prop_point
!
        IF ( nfun_reg(i,index) > 10 ) then
!
!                       General Code
!
           last = first + nfun_reg(i,index) - 1
           CALL v_so_mat_v_gen(v(first:last,1:nj),                               &
                               v_scr(first:last,1:nj),                            &
                               mat_reg(i,index)%exp_t_mat_z(:,:,prop_point))
        ELSE IF (nfun_reg(i,index) == 2) then
! 
!                       Special case for Two by Two
!
!           write(iout,*) 'Input Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
           last = first + 1
           CALL v_so_mat_v_2(v(first:last,1:nj),                                &
                             v_scr(first:last,1:nj),                            &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point))
!           write(iout,*) 'Output Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
        ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Special case for Three by Three
!
!           write(iout,*) 'Input Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
           last = first + 2
           CALL v_so_mat_v_3(v(first:last,1:nj),                               &
                             v_scr(first:last,1:nj),                           &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point))
!           write(iout,*) 'Output Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
        ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Special case for Four by Four
!
           last = first + 3              
           CALL v_so_mat_v_4(v(first:last,1:nj),                                &
                             v_scr(first:last,1:nj),                            &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point))
        ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Special case for Five by Five
!
           last = first + 4
           CALL v_so_mat_v_5(v(first:last,1:nj),                                &
                             v_scr(first:last,1:nj),                            &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point))
        ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Special case for Six by Six
!
!           write(iout,*) 'Input Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
           last = first + 5
           CALL v_so_mat_v_6(v(first:last,1:nj),                                &
                             v_scr(first:last,1:nj),                            &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point))
!           write(iout,*) 'Output Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
        ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Special case for Seven by Seven
!
           last = first + 6
           CALL v_so_mat_v_7(v(first:last,1:nj),                                &
                             v_scr(first:last,1:nj),                            &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point))
        ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Special case for Eight by Eight
!
           last = first + 7
           CALL v_so_mat_v_8(v(first:last,1:nj),                                &
                             v_scr(first:last,1:nj),                            &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point))
        ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Special case for Nine by Nine
!
           last = first + 8
           CALL v_so_mat_v_9(v(first:last,1:nj),                                &
                             v_scr(first:last,1:nj),                            &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point))
        ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Special case for Ten by Ten
!
           last = first + 9
           CALL v_so_mat_v_10(v(first:last,1:nj),                               &
                              v_scr(first:last,1:nj),                           &
                              mat_reg(i,index)%exp_t_mat_z(:,:,prop_point))
        END IF
        IF (i /= num_reg(index) ) then
            first = first + nfun_reg(i,index)                                   &
                          +                                                     &
                            nfun_reg(i+1,index)                                 &
                          -                                                     &
                            2
        END IF
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE exp_off_diag_m_v_z
!
!*deck exp_off_diag_v_m_2_z
!***begin prologue     exp_off_diag_v_m_2_z     
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
!***end prologue       exp_off_diag_v_m_2_z
!
  SUBROUTINE exp_off_diag_v_m_2_z(v,                                         &
                                  v_scr,                                     &
                                  ni,nj,index)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, index, ntrips
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  INTEGER                                  :: i, trips
  INTEGER                                  :: first, last
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
     first = locate(trips)
     DO i=begin(trips), num_reg(index), 2
        IF ( nfun_reg(i,index) > 10 ) then
!
!                        General Case
!
           last = first + nfun_reg(i,index) - 1
           CALL v_v_so_mat_gen(v(1:ni,first:last),                             &
                               v_scr(1:ni,first:last),                         &
                               mat_reg(i,index)%exp_t_mat_z(:,:,prop_point)) 
!
        ELSE IF (nfun_reg(i,index) == 2) then
!
!                       Two by Two
! 
           last= first + 1
           CALL v_v_so_mat_2(v(1:ni,first:last),                              &
                             v_scr(:,locate(trips):),                         &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point)) 
       ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Three by Three
!
           last= first + 2
           CALL v_v_so_mat_3(v(1:ni,first:last),                              &
                             v_scr(1:ni,first:last),                          &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point)) 
        ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Four by Four
!
           last= first + 3
           CALL v_v_so_mat_4(v(1:ni,first:last),                              &
                             v_scr(1:ni,first:last),                          &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point)) 
        ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Five by Five
!
           last= first + 4
           CALL v_v_so_mat_5(v(1:ni,first:last),                              &
                             v_scr(1:ni,first:last),                          &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point)) 
        ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Six by Six
!
           last= first + 5
           CALL v_v_so_mat_6(v(1:ni,first:last),                              &
                             v_scr(1:ni,first:last),                          &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point)) 
        ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Seven by Seven
!
           last= first + 6
           CALL v_v_so_mat_7(v(1:ni,first:last),                               &
                             v_scr(1:ni,first:last),                           &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point)) 
        ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Eight by Eight
!
           last= first + 7
           CALL v_v_so_mat_8(v(1:ni,first:last),                               &
                             v_scr(1:ni,first:last),                           &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point)) 
        ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Nine by Nine
!
           last= first + 8
           CALL v_v_so_mat_9(v(1:ni,first:last),                               &
                             v_scr(1:ni,first:last),                           &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point)) 
        ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Ten by Ten
!
           last= first + 9
           CALL v_v_so_mat_10(v(1:ni,first:last),                              &
                              v_scr(1:ni,first:last),                          &
                              mat_reg(i,index)%exp_t_mat_z(:,:,prop_point)) 
        END IF
        IF (i /= num_reg(index) ) then
            locate(trips) = locate(trips) + nfun_reg(i,index)                     &
                                          +                                       &
                            nfun_reg(i+1,index)                                   &
                                          -                                       &
                                          2
        END IF
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE exp_off_diag_v_m_2_z
!*deck exp_off_diag_v_m_3_z
!***begin prologue     exp_off_diag_v_m_3_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine simply calls
!***                   exp_off_diag_v_m_2 in a loop over a dummy
!***                   index, k
!
!***references
!***routines called
!***end prologue       exp_off_diag_v_m_3_z
!
  SUBROUTINE exp_off_diag_v_m_3_z(v,                                         &
                                  v_scr,                                     &
                                  ni,nj,nk,index)
  IMPLICIT NONE
  INTEGER                             :: ni, nj, nk, index
  COMPLEX*16, DIMENSION(ni,nj,nk)     :: v
  COMPLEX*16, DIMENSION(ni,nj,nk)     :: v_scr
  INTEGER                             :: k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  DO k=1,nk
     CALL exp_off_diag_v_m_2_z(v(:,:,k),                                     &
                               v_scr(:,:,k),                                 &
                               ni,nj,index) 
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE exp_off_diag_v_m_3_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END MODULE  split_operator_matrix_vector_multiply_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
