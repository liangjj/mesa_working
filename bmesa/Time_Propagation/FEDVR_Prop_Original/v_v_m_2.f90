!deck v_v_m_2
!***begin prologue     v_v_m_2     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagation matrix multiplies
!**                    for the special case of 2*2 matrices.
!
!***references
!***routines called
!***end prologue       v_v_m_2
!
  SUBROUTINE v_v_m_2(real_v,imag_v,real_v_scr,imag_v_scr, &
                     cosine_t_mat,sine_t_mat,ni,nj)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: real_v, imag_v
  REAL*8, DIMENSION(ni,nj)                 :: real_v_scr, imag_v_scr
  REAL*8, DIMENSION(2,2)                   :: cosine_t_mat
  REAL*8, DIMENSION(2,2)                   :: sine_t_mat
  INTEGER                                  :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  DO i=1,ni
     real_v_scr(i,1) = real_v(i,1) * cosine_t_mat(1,1)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,1)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,1)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,1)    
     real_v_scr(i,2) = real_v(i,1) * cosine_t_mat(1,2)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,2)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,2)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,2)    
     imag_v_scr(i,1) = imag_v(i,1) * cosine_t_mat(1,1)     &
                                      +                    &
                       imag_v(i,2) * cosine_t_mat(2,1)     &
                                      +                    &
                       real_v(i,1) * sine_t_mat(1,1)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,1) 
     imag_v_scr(i,2) = imag_v(i,1) * cosine_t_mat(1,2)     &
                                      +                    &
                       imag_v(i,2) * cosine_t_mat(2,2)     &
                                      +                    &
                       real_v(i,1) * sine_t_mat(1,2)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,2) 
  END DO
!
!       Copy the temporary vector back to the input vector.
!
  real_v(:,1:2) = real_v_scr(:,1:2) 
  imag_v(:,1:2) = imag_v_scr(:,1:2) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_m_2
