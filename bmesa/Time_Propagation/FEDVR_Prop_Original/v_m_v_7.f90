!deck v_m_v_7
!***begin prologue     v_m_v_7     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multipleis
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_m_v_7
!
  SUBROUTINE v_m_v_7(real_v,imag_v,real_v_scr,imag_v_scr, &
                     cosine_t_mat,sine_t_mat,ni,nj)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: real_v, imag_v
  REAL*8, DIMENSION(ni,nj)                 :: real_v_scr, imag_v_scr
  REAL*8, DIMENSION(7,7)                   :: cosine_t_mat
  REAL*8, DIMENSION(7,7)                   :: sine_t_mat
  INTEGER                                  :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  DO i=1,nj
     real_v_scr(1,i) = cosine_t_mat(1,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(1,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(1,3) * real_v(3,i)   &
                                         +               &   
                       cosine_t_mat(1,4) * real_v(4,i)   &
                                         +               &
                       cosine_t_mat(1,5) * real_v(5,i)   &   
                                         +               &
                       cosine_t_mat(1,6) * real_v(6,i)   &   
                                         +               &
                       cosine_t_mat(1,7) * real_v(7,i)   &   
                                         -               &
                       sine_t_mat(1,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(1,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(1,3)   * imag_v(3,i)   &
                                         -               &
                       sine_t_mat(1,4)   * imag_v(4,i)   &
                                         -               &
                       sine_t_mat(1,5)   * imag_v(5,i)   &
                                         -               &
                       sine_t_mat(1,6)   * imag_v(6,i)   &
                                         -               &
                       sine_t_mat(1,7)   * imag_v(7,i)   
     real_v_scr(2,i) = cosine_t_mat(2,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(2,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(2,3) * real_v(3,i)   &
                                         +               &   
                       cosine_t_mat(2,4) * real_v(4,i)   &
                                         +               &
                       cosine_t_mat(2,5) * real_v(5,i)   &   
                                         +               &
                       cosine_t_mat(2,6) * real_v(6,i)   &   
                                         +               &
                       cosine_t_mat(2,7) * real_v(7,i)   &   
                                         -               &
                       sine_t_mat(2,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(2,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(2,3)   * imag_v(3,i)   &
                                         -               &
                       sine_t_mat(2,4)   * imag_v(4,i)   &
                                         -               &
                       sine_t_mat(2,5)   * imag_v(5,i)   &
                                         -               &
                       sine_t_mat(2,6)   * imag_v(6,i)   &
                                         -               &
                       sine_t_mat(2,7)   * imag_v(7,i)   
     real_v_scr(3,i) = cosine_t_mat(3,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(3,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(3,3) * real_v(3,i)   &
                                         +               &   
                       cosine_t_mat(3,4) * real_v(4,i)   &
                                         +               &
                       cosine_t_mat(3,5) * real_v(5,i)   &   
                                         +               &
                       cosine_t_mat(3,6) * real_v(6,i)   &   
                                         +               &
                       cosine_t_mat(3,7) * real_v(7,i)   &   
                                         -               &
                       sine_t_mat(3,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(3,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(3,3)   * imag_v(3,i)   &
                                         -               &
                       sine_t_mat(3,4)   * imag_v(4,i)   &
                                         -               &
                       sine_t_mat(3,5)   * imag_v(5,i)   &
                                         -               &
                       sine_t_mat(3,6)   * imag_v(6,i)   &
                                         -               &
                       sine_t_mat(3,7)   * imag_v(7,i)   
     real_v_scr(4,i) = cosine_t_mat(4,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(4,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(4,3) * real_v(3,i)   &
                                         +               &   
                       cosine_t_mat(4,4) * real_v(4,i)   &
                                         +               &
                       cosine_t_mat(4,5) * real_v(5,i)   &   
                                         +               &
                       cosine_t_mat(4,6) * real_v(6,i)   &   
                                         +               &
                       cosine_t_mat(4,7) * real_v(7,i)   &   
                                         -               &
                       sine_t_mat(4,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(4,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(4,3)   * imag_v(3,i)   &
                                         -               &
                       sine_t_mat(4,4)   * imag_v(4,i)   &
                                         -               &
                       sine_t_mat(4,5)   * imag_v(5,i)   &
                                         -               &
                       sine_t_mat(4,6)   * imag_v(6,i)   &
                                         -               &
                       sine_t_mat(4,7)   * imag_v(7,i)   
     real_v_scr(5,i) = cosine_t_mat(5,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(5,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(5,3) * real_v(3,i)   &
                                         +               &   
                       cosine_t_mat(5,4) * real_v(4,i)   &
                                         +               &
                       cosine_t_mat(5,5) * real_v(5,i)   &   
                                         +               &
                       cosine_t_mat(5,6) * real_v(6,i)   &   
                                         +               &
                       cosine_t_mat(5,7) * real_v(7,i)   &   
                                         -               &
                       sine_t_mat(5,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(5,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(5,3)   * imag_v(3,i)   &
                                         -               &
                       sine_t_mat(5,4)   * imag_v(4,i)   &
                                         -               &
                       sine_t_mat(5,5)   * imag_v(5,i)   &
                                         -               &
                       sine_t_mat(5,6)   * imag_v(6,i)   &
                                         -               &
                       sine_t_mat(5,7)   * imag_v(7,i)   
     real_v_scr(6,i) = cosine_t_mat(6,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(6,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(6,3) * real_v(3,i)   &
                                         +               &   
                       cosine_t_mat(6,4) * real_v(4,i)   &
                                         +               &
                       cosine_t_mat(6,5) * real_v(5,i)   &   
                                         +               &
                       cosine_t_mat(6,6) * real_v(6,i)   &   
                                         +               &
                       cosine_t_mat(6,7) * real_v(7,i)   &   
                                         -               &
                       sine_t_mat(6,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(6,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(6,3)   * imag_v(3,i)   &
                                         -               &
                       sine_t_mat(6,4)   * imag_v(4,i)   &
                                         -               &
                       sine_t_mat(6,5)   * imag_v(5,i)   &
                                         -               &
                       sine_t_mat(6,6)   * imag_v(6,i)   &
                                         -               &
                       sine_t_mat(6,7)   * imag_v(7,i)   
     real_v_scr(7,i) = cosine_t_mat(7,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(7,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(7,3) * real_v(3,i)   &
                                         +               &   
                       cosine_t_mat(7,4) * real_v(4,i)   &
                                         +               &
                       cosine_t_mat(7,5) * real_v(5,i)   &   
                                         +               &
                       cosine_t_mat(7,6) * real_v(6,i)   &   
                                         +               &
                       cosine_t_mat(7,7) * real_v(7,i)   &   
                                         -               &
                       sine_t_mat(7,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(7,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(7,3)   * imag_v(3,i)   &
                                         -               &
                       sine_t_mat(7,4)   * imag_v(4,i)   &
                                         -               &
                       sine_t_mat(7,5)   * imag_v(5,i)   &
                                         -               &
                       sine_t_mat(7,6)   * imag_v(6,i)   &
                                         -               &
                       sine_t_mat(7,7)   * imag_v(7,i)   
     imag_v_scr(1,i) = cosine_t_mat(1,1) * imag_v(1,i)   &
                                         +               &
                       cosine_t_mat(1,2) * imag_v(2,i)   &
                                         +               &
                       cosine_t_mat(1,3) * imag_v(3,i)   &
                                         +               &
                       cosine_t_mat(1,4) * imag_v(4,i)   &
                                         +               &
                       cosine_t_mat(1,5) * imag_v(5,i)   &
                                         +               &
                       cosine_t_mat(1,6) * imag_v(6,i)   &
                                         +               &
                       cosine_t_mat(1,7) * imag_v(7,i)   &
                                         +               &
                       sine_t_mat(1,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(1,2)   * real_v(2,i)   &
                                         +               & 
                       sine_t_mat(1,3)   * real_v(3,i)   &
                                         +               &
                       sine_t_mat(1,4)   * real_v(4,i)   &
                                         +               &
                       sine_t_mat(1,5)   * real_v(5,i)   &
                                         +               &
                       sine_t_mat(1,6)   * real_v(6,i)   &
                                         +               &
                       sine_t_mat(1,7)   * real_v(7,i)   
     imag_v_scr(2,i) = cosine_t_mat(2,1) * imag_v(1,i)   &
                                         +               &
                       cosine_t_mat(2,2) * imag_v(2,i)   &
                                         +               &
                       cosine_t_mat(2,3) * imag_v(3,i)   &
                                         +               &
                       cosine_t_mat(2,4) * imag_v(4,i)   &
                                         +               &
                       cosine_t_mat(2,5) * imag_v(5,i)   &
                                         +               &
                       cosine_t_mat(2,6) * imag_v(6,i)   &
                                         +               &
                       cosine_t_mat(2,7) * imag_v(7,i)   &
                                         +               &
                       sine_t_mat(2,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(2,2)   * real_v(2,i)   &
                                         +               & 
                       sine_t_mat(2,3)   * real_v(3,i)   &
                                         +               &
                       sine_t_mat(2,4)   * real_v(4,i)   &
                                         +               &
                       sine_t_mat(2,5)   * real_v(5,i)   &
                                         +               &
                       sine_t_mat(2,6)   * real_v(6,i)   &
                                         +               &
                       sine_t_mat(2,7)   * real_v(7,i)   
     imag_v_scr(3,i) = cosine_t_mat(3,1) * imag_v(1,i)   &
                                         +               &
                       cosine_t_mat(3,2) * imag_v(2,i)   &
                                         +               &
                       cosine_t_mat(3,3) * imag_v(3,i)   &
                                         +               &
                       cosine_t_mat(3,4) * imag_v(4,i)   &
                                         +               &
                       cosine_t_mat(3,5) * imag_v(5,i)   &
                                         +               &
                       cosine_t_mat(3,6) * imag_v(6,i)   &
                                         +               &
                       cosine_t_mat(3,7) * imag_v(7,i)   &
                                         +               &
                       sine_t_mat(3,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(3,2)   * real_v(2,i)   &
                                         +               & 
                       sine_t_mat(3,3)   * real_v(3,i)   &
                                         +               &
                       sine_t_mat(3,4)   * real_v(4,i)   &
                                         +               &
                       sine_t_mat(3,5)   * real_v(5,i)   &
                                         +               &
                       sine_t_mat(3,6)   * real_v(6,i)   &
                                         +               &
                       sine_t_mat(3,7)   * real_v(7,i)   
     imag_v_scr(4,i) = cosine_t_mat(4,1) * imag_v(1,i)   &
                                         +               &
                       cosine_t_mat(4,2) * imag_v(2,i)   &
                                         +               &
                       cosine_t_mat(4,3) * imag_v(3,i)   &
                                         +               &
                       cosine_t_mat(4,4) * imag_v(4,i)   &
                                         +               &
                       cosine_t_mat(4,5) * imag_v(5,i)   &
                                         +               &
                       cosine_t_mat(4,6) * imag_v(6,i)   &
                                         +               &
                       cosine_t_mat(4,7) * imag_v(7,i)   &
                                         +               &
                       sine_t_mat(4,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(4,2)   * real_v(2,i)   &
                                         +               & 
                       sine_t_mat(4,3)   * real_v(3,i)   &
                                         +               &
                       sine_t_mat(4,4)   * real_v(4,i)   &
                                         +               &
                       sine_t_mat(4,5)   * real_v(5,i)   &
                                         +               &
                       sine_t_mat(4,6)   * real_v(6,i)   &
                                         +               &
                       sine_t_mat(4,7)   * real_v(7,i)   
     imag_v_scr(5,i) = cosine_t_mat(5,1) * imag_v(1,i)   &
                                         +               &
                       cosine_t_mat(5,2) * imag_v(2,i)   &
                                         +               &
                       cosine_t_mat(5,3) * imag_v(3,i)   &
                                         +               &
                       cosine_t_mat(5,4) * imag_v(4,i)   &
                                         +               &
                       cosine_t_mat(5,5) * imag_v(5,i)   &
                                         +               &
                       cosine_t_mat(5,6) * imag_v(6,i)   &
                                         +               &
                       cosine_t_mat(5,7) * imag_v(7,i)   &
                                         +               &
                       sine_t_mat(5,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(5,2)   * real_v(2,i)   &
                                         +               & 
                       sine_t_mat(5,3)   * real_v(3,i)   &
                                         +               &
                       sine_t_mat(5,4)   * real_v(4,i)   &
                                         +               &
                       sine_t_mat(5,5)   * real_v(5,i)   &
                                         +               &
                       sine_t_mat(5,6)   * real_v(6,i)   &
                                         +               &
                       sine_t_mat(5,7)   * real_v(7,i)   
     imag_v_scr(6,i) = cosine_t_mat(6,1) * imag_v(1,i)   &
                                         +               &
                       cosine_t_mat(6,2) * imag_v(2,i)   &
                                         +               &
                       cosine_t_mat(6,3) * imag_v(3,i)   &
                                         +               &
                       cosine_t_mat(6,4) * imag_v(4,i)   &
                                         +               &
                       cosine_t_mat(6,5) * imag_v(5,i)   &
                                         +               &
                       cosine_t_mat(6,6) * imag_v(6,i)   &
                                         +               &
                       cosine_t_mat(6,7) * imag_v(7,i)   &
                                         +               &
                       sine_t_mat(6,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(6,2)   * real_v(2,i)   &
                                         +               & 
                       sine_t_mat(6,3)   * real_v(3,i)   &
                                         +               &
                       sine_t_mat(6,4)   * real_v(4,i)   &
                                         +               &
                       sine_t_mat(6,5)   * real_v(5,i)   &
                                         +               &
                       sine_t_mat(6,6)   * real_v(6,i)   &
                                         +               &
                       sine_t_mat(6,7)   * real_v(7,i)   
     imag_v_scr(7,i) = cosine_t_mat(7,1) * imag_v(1,i)   &
                                         +               &
                       cosine_t_mat(7,2) * imag_v(2,i)   &
                                         +               &
                       cosine_t_mat(7,3) * imag_v(3,i)   &
                                         +               &
                       cosine_t_mat(7,4) * imag_v(4,i)   &
                                         +               &
                       cosine_t_mat(7,5) * imag_v(5,i)   &
                                         +               &
                       cosine_t_mat(7,6) * imag_v(6,i)   &
                                         +               &
                       cosine_t_mat(7,7) * imag_v(7,i)   &
                                         +               &
                       sine_t_mat(7,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(7,2)   * real_v(2,i)   &
                                         +               & 
                       sine_t_mat(7,3)   * real_v(3,i)   &
                                         +               &
                       sine_t_mat(7,4)   * real_v(4,i)   &
                                         +               &
                       sine_t_mat(7,5)   * real_v(5,i)   &
                                         +               &
                       sine_t_mat(7,6)   * real_v(6,i)   &
                                         +               &
                       sine_t_mat(7,7)   * real_v(7,i)   
  END DO
  real_v(1:7,:) = real_v_scr(1:7,:)
  imag_v(1:7,:) = imag_v_scr(1:7,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_m_v_7
