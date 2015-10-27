!deck v_m_v_2
!***begin prologue     v_m_v_2     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies for
!***                   the special case of 2*2 matrices.
!
!***references
!***routines called
!***end prologue       v_m_v_2
!
  SUBROUTINE v_m_v_2(real_v,imag_v,real_v_scr,imag_v_scr,    &
                     cosine_t_mat,sine_t_mat,ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: real_v, imag_v
  REAL*8, DIMENSION(ni,nj)                 :: real_v_scr, imag_v_scr
  REAL*8, DIMENSION(2,2)                   :: cosine_t_mat
  REAL*8, DIMENSION(2,2)                   :: sine_t_mat
  INTEGER                                  :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  DO i=1,nj
     real_v_scr(1,i) = cosine_t_mat(1,1) * real_v(1,i)     &
                                         +                 &
                       cosine_t_mat(1,2) * real_v(2,i)     &
                                         -                 &
                       sine_t_mat(1,1)   * imag_v(1,i)     &
                                         -                 &
                       sine_t_mat(1,2)   * imag_v(2,i)  
     real_v_scr(2,i) = cosine_t_mat(2,1) * real_v(1,i)     &
                                         +                 &
                       cosine_t_mat(2,2) * real_v(2,i)     &
                                         -                 &
                       sine_t_mat(2,1)   * imag_v(1,i)     &
                                         -                 &
                       sine_t_mat(2,2)   * imag_v(2,i)  
     imag_v_scr(1,i) = cosine_t_mat(1,1) * imag_v(1,i)     &
                                         +                 &
                       cosine_t_mat(1,2) * imag_v(2,i)     &
                                         +                 &
                       sine_t_mat(1,1)   * real_v(1,i)     &
                                         +                 &
                       sine_t_mat(1,2)   * real_v(2,i)  
     imag_v_scr(2,i) = cosine_t_mat(2,1) * imag_v(1,i)     &
                                         +                 &
                       cosine_t_mat(2,2) * imag_v(2,i)     &
                                         +                 &
                       sine_t_mat(2,1)   * real_v(1,i)     &
                                         +                 &
                       sine_t_mat(2,2)   * real_v(2,i)  
  END DO
  real_v(1:2,:) = real_v_scr(1:2,:)
  imag_v(1:2,:) = imag_v_scr(1:2,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_m_v_2
