!deck v_m_v_3
!***begin prologue     v_m_v_3     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multipleis
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_m_v_3
!
  SUBROUTINE v_m_v_3(real_v,imag_v,real_v_scr,imag_v_scr, &
                     cosine_t_mat,sine_t_mat,ni,nj)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: real_v, imag_v
  REAL*8, DIMENSION(ni,nj)                 :: real_v_scr, imag_v_scr
  REAL*8, DIMENSION(3,3)                   :: cosine_t_mat
  REAL*8, DIMENSION(3,3)                   :: sine_t_mat
  INTEGER                                  :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  DO i=1,nj
     real_v_scr(1,i) = cosine_t_mat(1,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(1,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(1,3) * real_v(3,i)   &
                                         -               &
                       sine_t_mat(1,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(1,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(1,3)   * imag_v(3,i)  
     real_v_scr(2,i) = cosine_t_mat(2,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(2,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(2,3) * real_v(3,i)   &
                                         -               &
                       sine_t_mat(2,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(2,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(2,3)   * imag_v(3,i)  
     real_v_scr(3,i) = cosine_t_mat(3,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(3,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(3,3) * real_v(3,i)   &
                                         -               &
                       sine_t_mat(3,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(3,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(3,3)   * imag_v(3,i)  
     imag_v_scr(1,i) = cosine_t_mat(1,1) * imag_v(1,i)   &
                                         +               &
                       cosine_t_mat(1,2) * imag_v(2,i)   &
                                         +               &
                       cosine_t_mat(1,3) * imag_v(3,i)   &
                                         +               &
                       sine_t_mat(1,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(1,2)   * real_v(2,i)   &
                                         +               & 
                       sine_t_mat(1,3)   * real_v(3,i)  
     imag_v_scr(2,i) = cosine_t_mat(2,1) * imag_v(1,i)   &
                                         +               &
                       cosine_t_mat(2,2) * imag_v(2,i)   &
                                         +               &
                       cosine_t_mat(2,3) * imag_v(3,i)   &
                                         +               &
                       sine_t_mat(2,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(2,2)   * real_v(2,i)   &
                                         +               &
                       sine_t_mat(2,3)   * real_v(3,i)  
     imag_v_scr(3,i) = cosine_t_mat(3,1) * imag_v(1,i)   &
                                         +               &
                       cosine_t_mat(3,2) * imag_v(2,i)   &
                                         +               &
                       cosine_t_mat(3,3) * imag_v(3,i)   &
                                         +               &
                       sine_t_mat(3,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(3,2)   * real_v(2,i)   &
                                         +               &
                       sine_t_mat(3,3)   * real_v(3,i)  
  END DO
  real_v(1:3,:) = real_v_scr(1:3,:)
  imag_v(1:3,:) = imag_v_scr(1:3,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_m_v_3
