!deck v_m_v_9
!***begin prologue     v_m_v_9     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 6*6 matrices.
!
!***references
!***routines called
!***end prologue       v_m_v_9
!
  SUBROUTINE v_m_v_9(real_v,imag_v,real_v_scr,imag_v_scr, &
                     cosine_t_mat,sine_t_mat,ni,nj)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: real_v, imag_v
  REAL*8, DIMENSION(ni,nj)                 :: real_v_scr, imag_v_scr
  REAL*8, DIMENSION(9,9)                   :: cosine_t_mat
  REAL*8, DIMENSION(9,9)                   :: sine_t_mat
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
                                         +               &
                       cosine_t_mat(1,8) * real_v(8,i)   &   
                                         +               &
                       cosine_t_mat(1,9) * real_v(9,i)   &   
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
                       sine_t_mat(1,7)   * imag_v(7,i)   &
                                         -               &
                       sine_t_mat(1,8)   * imag_v(8,i)   &
                                         -               &
                       sine_t_mat(1,9)   * imag_v(9,i)   
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
                                         +               &
                       cosine_t_mat(2,8) * real_v(8,i)   &   
                                         +               &
                       cosine_t_mat(2,9) * real_v(9,i)   &   
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
                       sine_t_mat(2,7)   * imag_v(7,i)   &
                                         -               &
                       sine_t_mat(2,8)   * imag_v(8,i)   &
                                         -               &
                       sine_t_mat(2,9)   * imag_v(9,i)   
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
                                         +               &
                       cosine_t_mat(3,8) * real_v(8,i)   &   
                                         +               &
                       cosine_t_mat(3,9) * real_v(9,i)   &   
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
                       sine_t_mat(3,7)   * imag_v(7,i)   &
                                         -               &
                       sine_t_mat(3,8)   * imag_v(8,i)   &
                                         -               &
                       sine_t_mat(3,9)   * imag_v(9,i)   
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
                                         +               &
                       cosine_t_mat(4,8) * real_v(8,i)   &   
                                         +               &
                       cosine_t_mat(4,9) * real_v(9,i)   &   
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
                       sine_t_mat(4,7)   * imag_v(7,i)   &
                                         -               &
                       sine_t_mat(4,8)   * imag_v(8,i)   &
                                         -               &
                       sine_t_mat(4,9)   * imag_v(9,i)   
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
                                         +               &
                       cosine_t_mat(5,8) * real_v(8,i)   &   
                                         +               &
                       cosine_t_mat(5,9) * real_v(9,i)   &   
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
                       sine_t_mat(5,7)   * imag_v(7,i)   &
                                         -               &
                       sine_t_mat(5,8)   * imag_v(8,i)   &
                                         -               &
                       sine_t_mat(5,9)   * imag_v(9,i)   
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
                                         +               &
                       cosine_t_mat(6,8) * real_v(8,i)   &   
                                         +               &
                       cosine_t_mat(6,9) * real_v(9,i)   &   
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
                       sine_t_mat(6,7)   * imag_v(7,i)   &
                                         -               &
                       sine_t_mat(6,8)   * imag_v(8,i)   &
                                         -               &
                       sine_t_mat(6,9)   * imag_v(9,i)   
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
                                         +               &
                       cosine_t_mat(7,8) * real_v(8,i)   &   
                                         +               &
                       cosine_t_mat(7,9) * real_v(9,i)   &   
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
                       sine_t_mat(7,7)   * imag_v(7,i)   &
                                         -               &
                       sine_t_mat(7,8)   * imag_v(8,i)   &
                                         -               &
                       sine_t_mat(7,9)   * imag_v(9,i)   
     real_v_scr(8,i) = cosine_t_mat(8,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(8,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(8,3) * real_v(3,i)   &
                                         +               &   
                       cosine_t_mat(8,4) * real_v(4,i)   &
                                         +               &
                       cosine_t_mat(8,5) * real_v(5,i)   &   
                                         +               &
                       cosine_t_mat(8,6) * real_v(6,i)   &   
                                         +               &
                       cosine_t_mat(8,7) * real_v(7,i)   &   
                                         +               &
                       cosine_t_mat(8,8) * real_v(8,i)   &   
                                         +               &
                       cosine_t_mat(8,9) * real_v(9,i)   &   
                                         -               &
                       sine_t_mat(8,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(8,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(8,3)   * imag_v(3,i)   &
                                         -               &
                       sine_t_mat(8,4)   * imag_v(4,i)   &
                                         -               &
                       sine_t_mat(8,5)   * imag_v(5,i)   &
                                         -               &
                       sine_t_mat(8,6)   * imag_v(6,i)   &
                                         -               &
                       sine_t_mat(8,7)   * imag_v(7,i)   &
                                         -               &
                       sine_t_mat(8,8)   * imag_v(8,i)   &
                                         -               &
                       sine_t_mat(8,9)   * imag_v(9,i)   
     real_v_scr(9,i) = cosine_t_mat(9,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(9,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(9,3) * real_v(3,i)   &
                                         +               &   
                       cosine_t_mat(9,4) * real_v(4,i)   &
                                         +               &
                       cosine_t_mat(9,5) * real_v(5,i)   &   
                                         +               &
                       cosine_t_mat(9,6) * real_v(6,i)   &   
                                         +               &
                       cosine_t_mat(9,7) * real_v(7,i)   &   
                                         +               &
                       cosine_t_mat(9,8) * real_v(8,i)   &   
                                         +               &
                       cosine_t_mat(9,9) * real_v(9,i)   &   
                                         -               &
                       sine_t_mat(9,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(9,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(9,3)   * imag_v(3,i)   &
                                         -               &
                       sine_t_mat(9,4)   * imag_v(4,i)   &
                                         -               &
                       sine_t_mat(9,5)   * imag_v(5,i)   &
                                         -               &
                       sine_t_mat(9,6)   * imag_v(6,i)   &
                                         -               &
                       sine_t_mat(9,7)   * imag_v(7,i)   &
                                         -               &
                       sine_t_mat(9,8)   * imag_v(8,i)   &
                                         -               &
                       sine_t_mat(9,9)   * imag_v(9,i)   

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
                       cosine_t_mat(1,8) * imag_v(8,i)   &
                                         +               &
                       cosine_t_mat(1,9) * imag_v(9,i)   &
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
                       sine_t_mat(1,7)   * real_v(7,i)   &
                                         +               &
                       sine_t_mat(1,8)   * real_v(8,i)   &
                                         +               &
                       sine_t_mat(1,9)   * real_v(9,i)   
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
                       cosine_t_mat(2,8) * imag_v(8,i)   &
                                         +               &
                       cosine_t_mat(2,9) * imag_v(9,i)   &
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
                       sine_t_mat(2,7)   * real_v(7,i)   &
                                         +               &
                       sine_t_mat(2,8)   * real_v(8,i)   &
                                         +               &
                       sine_t_mat(2,9)   * real_v(9,i)   
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
                       cosine_t_mat(3,8) * imag_v(8,i)   &
                                         +               &
                       cosine_t_mat(3,9) * imag_v(9,i)   &
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
                       sine_t_mat(3,7)   * real_v(7,i)   &
                                         +               &
                       sine_t_mat(3,8)   * real_v(8,i)   &
                                         +               &
                       sine_t_mat(3,9)   * real_v(9,i)   
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
                       cosine_t_mat(4,8) * imag_v(8,i)   &
                                         +               &
                       cosine_t_mat(4,9) * imag_v(9,i)   &
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
                       sine_t_mat(4,7)   * real_v(7,i)   &
                                         +               &
                       sine_t_mat(4,8)   * real_v(8,i)   &
                                         +               &
                       sine_t_mat(4,9)   * real_v(9,i)   
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
                       cosine_t_mat(5,8) * imag_v(8,i)   &
                                         +               &
                       cosine_t_mat(5,9) * imag_v(9,i)   &
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
                       sine_t_mat(5,7)   * real_v(7,i)   &
                                         +               &
                       sine_t_mat(5,8)   * real_v(8,i)   &
                                         +               &
                       sine_t_mat(5,9)   * real_v(9,i)   
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
                       cosine_t_mat(6,8) * imag_v(8,i)   &
                                         +               &
                       cosine_t_mat(6,9) * imag_v(9,i)   &
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
                       sine_t_mat(6,7)   * real_v(7,i)   &
                                         +               &
                       sine_t_mat(6,8)   * real_v(8,i)   &
                                         +               &
                       sine_t_mat(6,9)   * real_v(9,i)   
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
                       cosine_t_mat(7,8) * imag_v(8,i)   &
                                         +               &
                       cosine_t_mat(7,9) * imag_v(9,i)   &
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
                       sine_t_mat(7,7)   * real_v(7,i)   &
                                         +               &
                       sine_t_mat(7,8)   * real_v(8,i)   &
                                         +               &
                       sine_t_mat(7,9)   * real_v(9,i)   
     imag_v_scr(8,i) = cosine_t_mat(8,1) * imag_v(1,i)   &
                                         +               &
                       cosine_t_mat(8,2) * imag_v(2,i)   &
                                         +               &
                       cosine_t_mat(8,3) * imag_v(3,i)   &
                                         +               &
                       cosine_t_mat(8,4) * imag_v(4,i)   &
                                         +               &
                       cosine_t_mat(8,5) * imag_v(5,i)   &
                                         +               &
                       cosine_t_mat(8,6) * imag_v(6,i)   &
                                         +               &
                       cosine_t_mat(8,7) * imag_v(7,i)   &
                                         +               &
                       cosine_t_mat(8,8) * imag_v(8,i)   &
                                         +               &
                       cosine_t_mat(8,9) * imag_v(9,i)   &
                                         +               &
                       sine_t_mat(8,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(8,2)   * real_v(2,i)   &
                                         +               & 
                       sine_t_mat(8,3)   * real_v(3,i)   &
                                         +               &
                       sine_t_mat(8,4)   * real_v(4,i)   &
                                         +               &
                       sine_t_mat(8,5)   * real_v(5,i)   &
                                         +               &
                       sine_t_mat(8,6)   * real_v(6,i)   &
                                         +               &
                       sine_t_mat(8,7)   * real_v(7,i)   &
                                         +               &
                       sine_t_mat(8,8)   * real_v(8,i)   &
                                         +               &
                       sine_t_mat(8,9)   * real_v(9,i)   
     imag_v_scr(9,i) = cosine_t_mat(9,1) * imag_v(1,i)   &
                                         +               &
                       cosine_t_mat(9,2) * imag_v(2,i)   &
                                         +               &
                       cosine_t_mat(9,3) * imag_v(3,i)   &
                                         +               &
                       cosine_t_mat(9,4) * imag_v(4,i)   &
                                         +               &
                       cosine_t_mat(9,5) * imag_v(5,i)   &
                                         +               &
                       cosine_t_mat(9,6) * imag_v(6,i)   &
                                         +               &
                       cosine_t_mat(9,7) * imag_v(7,i)   &
                                         +               &
                       cosine_t_mat(9,8) * imag_v(8,i)   &
                                         +               &
                       cosine_t_mat(9,9) * imag_v(9,i)   &
                                         +               &
                       sine_t_mat(9,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(9,2)   * real_v(2,i)   &
                                         +               & 
                       sine_t_mat(9,3)   * real_v(3,i)   &
                                         +               &
                       sine_t_mat(9,4)   * real_v(4,i)   &
                                         +               &
                       sine_t_mat(9,5)   * real_v(5,i)   &
                                         +               &
                       sine_t_mat(9,6)   * real_v(6,i)   &
                                         +               &
                       sine_t_mat(9,7)   * real_v(7,i)   &
                                         +               &
                       sine_t_mat(9,8)   * real_v(8,i)   &
                                         +               &
                       sine_t_mat(9,9)   * real_v(9,i)   

  END DO
  real_v(1:9,:) = real_v_scr(1:9,:)
  imag_v(1:9,:) = imag_v_scr(1:9,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_m_v_9
