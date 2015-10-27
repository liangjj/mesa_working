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
!                 Real(V) = cos_t_mat * Real(V) - sin_t_mat * Imag(V)
!                 Imag(V) = cos_t_mat * Imag(V) + sin_t_mat * Real(V)
!
!------------------------------------------------------------------------------------
!***                   where cos_t_mat and sine_t_mat are the regional
!***                   propagators for a FEDVR or FD Hamiltonian.
!
!***                   The parameter nj is a dummy and can take on values
!***                   consistent with any dimensional problem.  In a 1D
!***                   case nj=1, in 2D nj=nx and in 3D nj=ny*nx.

!***references
!***routines called    v_m_v_gen, v_m_v_2, v_m_v_3
!
!***end prologue       exp_off_diag_m_v                     
!
  SUBROUTINE exp_off_diag_m_v(real_v,           &
                              imag_v,           &
                              real_v_scr,       &
                              imag_v_scr,       &
                              ni,nj,index)
  USE io
  USE dvr_global,                  ONLY    : num_reg, nfun_reg
  USE dvrprop_global,              ONLY    : mat_reg_d
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, index
  REAL*8, DIMENSION(ni,nj)                 :: real_v, imag_v
  REAL*8, DIMENSION(ni,nj)                 :: real_v_scr, imag_v_scr
  INTEGER                                  :: i, trips
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
  DO trips=1,3
     DO i = begin(trips), num_reg(index), 2
!
        IF ( nfun_reg(i,index) > 7 ) then
!
!                       General Code
!
           CALL v_m_v_gen(real_v(locate(trips),1),              &
                          imag_v(locate(trips),1),              &
                          real_v_scr(locate(trips),1),          &
                          imag_v_scr(locate(trips),1),          &
                          mat_reg_d(i,index)%cosine_t_mat,      &
                          mat_reg_d(i,index)%sine_t_mat,        &
                          ni,nj,nfun_reg(i,index))
        ELSE IF (nfun_reg(i,index) == 2) then
! 
!                       Special case for Two by Two
!
           CALL v_m_v_2(real_v(locate(trips),1),                &
                        imag_v(locate(trips),1),                &
                        real_v_scr(locate(trips),1),            &
                        imag_v_scr(locate(trips),1),            &
                        mat_reg_d(i,index)%cosine_t_mat,        &
                        mat_reg_d(i,index)%sine_t_mat,          &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Special case for Three by Three
!
           CALL v_m_v_3(real_v(locate(trips),1),                &
                        imag_v(locate(trips),1),                &
                        real_v_scr(locate(trips),1),            &
                        imag_v_scr(locate(trips),1),            &
                        mat_reg_d(i,index)%cosine_t_mat,        &
                        mat_reg_d(i,index)%sine_t_mat,          &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Special case for Four by Four
!
           CALL v_m_v_4(real_v(locate(trips),1),                &
                        imag_v(locate(trips),1),                &
                        real_v_scr(locate(trips),1),            &
                        imag_v_scr(locate(trips),1),            &
                        mat_reg_d(i,index)%cosine_t_mat,        &
                        mat_reg_d(i,index)%sine_t_mat,          &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Special case for Five by Five
!
           CALL v_m_v_5(real_v(locate(trips),1),                &
                        imag_v(locate(trips),1),                &
                        real_v_scr(locate(trips),1),            &
                        imag_v_scr(locate(trips),1),            &
                        mat_reg_d(i,index)%cosine_t_mat,        &
                        mat_reg_d(i,index)%sine_t_mat,          &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Special case for Six by Six
!
           CALL v_m_v_6(real_v(locate(trips),1),                &
                        imag_v(locate(trips),1),                &
                        real_v_scr(locate(trips),1),            &
                        imag_v_scr(locate(trips),1),            &
                        mat_reg_d(i,index)%cosine_t_mat,        &
                        mat_reg_d(i,index)%sine_t_mat,          &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Special case for Seven by Seven
!
           CALL v_m_v_7(real_v(locate(trips),1),                &
                        imag_v(locate(trips),1),                &
                        real_v_scr(locate(trips),1),            &
                        imag_v_scr(locate(trips),1),            &
                        mat_reg_d(i,index)%cosine_t_mat,        &
                        mat_reg_d(i,index)%sine_t_mat,          &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Special case for Eight by Eight
!
           CALL v_m_v_8(real_v(locate(trips),1),                &
                        imag_v(locate(trips),1),                &
                        real_v_scr(locate(trips),1),            &
                        imag_v_scr(locate(trips),1),            &
                        mat_reg_d(i,index)%cosine_t_mat,        &
                        mat_reg_d(i,index)%sine_t_mat,          &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Special case for Nine by Nine
!
           CALL v_m_v_9(real_v(locate(trips),1),                &
                        imag_v(locate(trips),1),                &
                        real_v_scr(locate(trips),1),            &
                        imag_v_scr(locate(trips),1),            &
                        mat_reg_d(i,index)%cosine_t_mat,        &
                        mat_reg_d(i,index)%sine_t_mat,          &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Special case for Ten by Ten
!
           CALL v_m_v_10(real_v(locate(trips),1),               &
                         imag_v(locate(trips),1),               &
                         real_v_scr(locate(trips),1),           &
                         imag_v_scr(locate(trips),1),           &
                         mat_reg_d(i,index)%cosine_t_mat,       &
                         mat_reg_d(i,index)%sine_t_mat,         &
                         ni,nj)
        END IF
        IF (i /= num_reg(index) ) then
            locate(trips) = locate(trips) + nfun_reg(i,index)   &
                                          +                     &
                            nfun_reg(i+1,index)                 &
                                          -                     &
                                          2
        END IF
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE exp_off_diag_m_v
