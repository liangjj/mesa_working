!deck v_v_m_8
!***begin prologue     v_v_m_8     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_v_m_8
!
  SUBROUTINE v_v_m_8(real_v,imag_v,real_v_scr,imag_v_scr, &
                     cosine_t_mat,sine_t_mat,ni,nj)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: real_v, imag_v
  REAL*8, DIMENSION(ni,nj)                 :: real_v_scr, imag_v_scr
  REAL*8, DIMENSION(8,8)                   :: cosine_t_mat
  REAL*8, DIMENSION(8,8)                   :: sine_t_mat
  INTEGER                                  :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  DO i=1,ni
     real_v_scr(i,1) = real_v(i,1) * cosine_t_mat(1,1)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,1)     &
                                      +                    &
                       real_v(i,3) * cosine_t_mat(3,1)     &
                                      +                    &
                       real_v(i,4) * cosine_t_mat(4,1)     &
                                      +                    &
                       real_v(i,5) * cosine_t_mat(5,1)     &
                                      +                    &
                       real_v(i,6) * cosine_t_mat(6,1)     &
                                      +                    &
                       real_v(i,7) * cosine_t_mat(7,1)     &
                                      +                    &
                       real_v(i,8) * cosine_t_mat(8,1)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,1)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,1)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,1)       &
                                      -                    &
                       imag_v(i,4) * sine_t_mat(4,1)       &
                                      -                    &
                       imag_v(i,5) * sine_t_mat(5,1)       &
                                      -                    &
                       imag_v(i,6) * sine_t_mat(6,1)       &
                                      -                    &
                       imag_v(i,7) * sine_t_mat(7,1)       &
                                      -                    &
                       imag_v(i,8) * sine_t_mat(8,1)    
     real_v_scr(i,2) = real_v(i,1) * cosine_t_mat(1,2)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,2)     &
                                      +                    &
                       real_v(i,3) * cosine_t_mat(3,2)     &
                                      +                    &
                       real_v(i,4) * cosine_t_mat(4,2)     &
                                      +                    &
                       real_v(i,5) * cosine_t_mat(5,2)     &
                                      +                    &
                       real_v(i,6) * cosine_t_mat(6,2)     &
                                      +                    &
                       real_v(i,7) * cosine_t_mat(7,2)     &
                                      +                    &
                       real_v(i,8) * cosine_t_mat(8,2)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,2)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,2)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,2)       &
                                      -                    &
                       imag_v(i,4) * sine_t_mat(4,2)       &
                                      -                    &
                       imag_v(i,5) * sine_t_mat(5,2)       &
                                      -                    &
                       imag_v(i,6) * sine_t_mat(6,2)       &
                                      -                    &
                       imag_v(i,7) * sine_t_mat(7,2)       &
                                      -                    &
                       imag_v(i,8) * sine_t_mat(8,2)    
     real_v_scr(i,3) = real_v(i,1) * cosine_t_mat(1,3)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,3)     &
                                      +                    &
                       real_v(i,3) * cosine_t_mat(3,3)     &
                                      +                    &
                       real_v(i,4) * cosine_t_mat(4,3)     &
                                      +                    &
                       real_v(i,5) * cosine_t_mat(5,3)     &
                                      +                    &
                       real_v(i,6) * cosine_t_mat(6,3)     &
                                      +                    &
                       real_v(i,7) * cosine_t_mat(7,3)     &
                                      +                    &
                       real_v(i,8) * cosine_t_mat(8,3)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,3)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,3)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,3)       &
                                      -                    &
                       imag_v(i,4) * sine_t_mat(4,3)       &
                                      -                    &
                       imag_v(i,5) * sine_t_mat(5,3)       &
                                      -                    &
                       imag_v(i,6) * sine_t_mat(6,3)       &
                                      -                    &
                       imag_v(i,7) * sine_t_mat(7,3)       &
                                      -                    &
                       imag_v(i,8) * sine_t_mat(8,3)    
     real_v_scr(i,4) = real_v(i,1) * cosine_t_mat(1,4)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,4)     &
                                      +                    &
                       real_v(i,3) * cosine_t_mat(3,4)     &
                                      +                    &
                       real_v(i,4) * cosine_t_mat(4,4)     &
                                      +                    &
                       real_v(i,5) * cosine_t_mat(5,4)     &
                                      +                    &
                       real_v(i,6) * cosine_t_mat(6,4)     &
                                      +                    &
                       real_v(i,7) * cosine_t_mat(7,4)     &
                                      +                    &
                       real_v(i,8) * cosine_t_mat(8,4)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,4)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,4)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,4)       &
                                      -                    &
                       imag_v(i,4) * sine_t_mat(4,4)       &
                                      -                    &
                       imag_v(i,5) * sine_t_mat(5,4)       &
                                      -                    &
                       imag_v(i,6) * sine_t_mat(6,4)       &
                                      -                    &
                       imag_v(i,7) * sine_t_mat(7,4)       &
                                      -                    &
                       imag_v(i,8) * sine_t_mat(8,4)    
     real_v_scr(i,5) = real_v(i,1) * cosine_t_mat(1,5)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,5)     &
                                      +                    &
                       real_v(i,3) * cosine_t_mat(3,5)     &
                                      +                    &
                       real_v(i,4) * cosine_t_mat(4,5)     &
                                      +                    &
                       real_v(i,5) * cosine_t_mat(5,5)     &
                                      +                    &
                       real_v(i,6) * cosine_t_mat(6,5)     &
                                      +                    &
                       real_v(i,7) * cosine_t_mat(7,5)     &
                                      +                    &
                       real_v(i,8) * cosine_t_mat(8,5)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,5)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,5)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,5)       &
                                      -                    &
                       imag_v(i,4) * sine_t_mat(4,5)       &
                                      -                    &
                       imag_v(i,5) * sine_t_mat(5,5)       &
                                      -                    &
                       imag_v(i,6) * sine_t_mat(6,5)       &
                                      -                    &
                       imag_v(i,7) * sine_t_mat(7,5)       &
                                      -                    &
                       imag_v(i,8) * sine_t_mat(8,1)    
     real_v_scr(i,6) = real_v(i,1) * cosine_t_mat(1,6)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,6)     &
                                      +                    &
                       real_v(i,3) * cosine_t_mat(3,6)     &
                                      +                    &
                       real_v(i,4) * cosine_t_mat(4,6)     &
                                      +                    &
                       real_v(i,5) * cosine_t_mat(5,6)     &
                                      +                    &
                       real_v(i,6) * cosine_t_mat(6,6)     &
                                      +                    &
                       real_v(i,7) * cosine_t_mat(7,6)     &
                                      +                    &
                       real_v(i,8) * cosine_t_mat(8,6)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,6)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,6)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,6)       &
                                      -                    &
                       imag_v(i,4) * sine_t_mat(4,6)       &
                                      -                    &
                       imag_v(i,5) * sine_t_mat(5,6)       &
                                      -                    &
                       imag_v(i,6) * sine_t_mat(6,6)       &
                                      -                    &
                       imag_v(i,7) * sine_t_mat(7,6)       &
                                      -                    &
                       imag_v(i,8) * sine_t_mat(8,6)    
     real_v_scr(i,7) = real_v(i,1) * cosine_t_mat(1,7)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,7)     &
                                      +                    &
                       real_v(i,3) * cosine_t_mat(3,7)     &
                                      +                    &
                       real_v(i,4) * cosine_t_mat(4,7)     &
                                      +                    &
                       real_v(i,5) * cosine_t_mat(5,7)     &
                                      +                    &
                       real_v(i,6) * cosine_t_mat(6,7)     &
                                      +                    &
                       real_v(i,7) * cosine_t_mat(7,7)     &
                                      +                    &
                       real_v(i,8) * cosine_t_mat(8,7)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,7)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,7)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,7)       &
                                      -                    &
                       imag_v(i,4) * sine_t_mat(4,7)       &
                                      -                    &
                       imag_v(i,5) * sine_t_mat(5,7)       &
                                      -                    &
                       imag_v(i,6) * sine_t_mat(6,7)       &
                                      -                    &
                       imag_v(i,7) * sine_t_mat(7,7)       &
                                      -                    &
                       imag_v(i,8) * sine_t_mat(8,7)    
     real_v_scr(i,8) = real_v(i,1) * cosine_t_mat(1,8)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,8)     &
                                      +                    &
                       real_v(i,3) * cosine_t_mat(3,8)     &
                                      +                    &
                       real_v(i,4) * cosine_t_mat(4,8)     &
                                      +                    &
                       real_v(i,5) * cosine_t_mat(5,8)     &
                                      +                    &
                       real_v(i,6) * cosine_t_mat(6,8)     &
                                      +                    &
                       real_v(i,7) * cosine_t_mat(7,8)     &
                                      +                    &
                       real_v(i,8) * cosine_t_mat(8,8)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,8)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,8)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,8)       &
                                      -                    &
                       imag_v(i,4) * sine_t_mat(4,8)       &
                                      -                    &
                       imag_v(i,5) * sine_t_mat(5,8)       &
                                      -                    &
                       imag_v(i,6) * sine_t_mat(6,8)       &
                                      -                    &
                       imag_v(i,7) * sine_t_mat(7,8)       &
                                      -                    &
                       imag_v(i,8) * sine_t_mat(8,8)    
     imag_v_scr(i,1) = imag_v(i,1) * cosine_t_mat(1,1)     &
                                      +                    &
                       imag_v(i,2) * cosine_t_mat(2,1)     &
                                      +                    &
                       imag_v(i,3) * cosine_t_mat(3,1)     &
                                      +                    &
                       imag_v(i,4) * cosine_t_mat(4,1)     &
                                      +                    &
                       imag_v(i,5) * cosine_t_mat(5,1)     &
                                      +                    &
                       imag_v(i,6) * cosine_t_mat(6,1)     &
                                      +                    &
                       imag_v(i,7) * cosine_t_mat(7,1)     &
                                      +                    &
                       imag_v(i,8) * cosine_t_mat(8,1)     &
                                      +                    &
                       real_v(i,1) * sine_t_mat(1,1)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,1)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,1)       &
                                      +                    &
                       real_v(i,4) * sine_t_mat(4,1)       &
                                      +                    &
                       real_v(i,5) * sine_t_mat(5,1)       &
                                      +                    &
                       real_v(i,6) * sine_t_mat(6,1)       &
                                      +                    &
                       real_v(i,7) * sine_t_mat(7,1)       &
                                      +                    &
                       real_v(i,8) * sine_t_mat(8,1) 
     imag_v_scr(i,2) = imag_v(i,1) * cosine_t_mat(1,2)     &
                                      +                    &
                       imag_v(i,2) * cosine_t_mat(2,2)     &
                                      +                    &
                       imag_v(i,3) * cosine_t_mat(3,2)     &
                                      +                    &
                       imag_v(i,4) * cosine_t_mat(4,2)     &
                                      +                    &
                       imag_v(i,5) * cosine_t_mat(5,2)     &
                                      +                    &
                       imag_v(i,6) * cosine_t_mat(6,2)     &
                                      +                    &
                       imag_v(i,7) * cosine_t_mat(7,2)     &
                                      +                    &
                       imag_v(i,8) * cosine_t_mat(8,2)     &
                                      +                    &
                       real_v(i,1) * sine_t_mat(1,2)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,2)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,2)       &
                                      +                    &
                       real_v(i,4) * sine_t_mat(4,2)       &
                                      +                    &
                       real_v(i,5) * sine_t_mat(5,2)       &
                                      +                    &
                       real_v(i,6) * sine_t_mat(6,2)       &
                                      +                    &
                       real_v(i,7) * sine_t_mat(7,2)       &
                                      +                    &
                       real_v(i,8) * sine_t_mat(8,2) 
     imag_v_scr(i,3) = imag_v(i,1) * cosine_t_mat(1,3)     &
                                      +                    &
                       imag_v(i,2) * cosine_t_mat(2,3)     &
                                      +                    &
                       imag_v(i,3) * cosine_t_mat(3,3)     &
                                      +                    &
                       imag_v(i,4) * cosine_t_mat(4,3)     &
                                      +                    &
                       imag_v(i,5) * cosine_t_mat(5,3)     &
                                      +                    &
                       imag_v(i,6) * cosine_t_mat(6,3)     &
                                      +                    &
                       imag_v(i,7) * cosine_t_mat(7,3)     &
                                      +                    &
                       imag_v(i,8) * cosine_t_mat(8,3)     &
                                      +                    &
                       real_v(i,1) * sine_t_mat(1,3)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,3)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,3)       &
                                      +                    &
                       real_v(i,4) * sine_t_mat(4,3)       &
                                      +                    &
                       real_v(i,5) * sine_t_mat(5,3)       &
                                      +                    &
                       real_v(i,6) * sine_t_mat(6,3)       &
                                      +                    &
                       real_v(i,7) * sine_t_mat(7,3)       &
                                      +                    &
                       real_v(i,8) * sine_t_mat(8,3) 
     imag_v_scr(i,4) = imag_v(i,1) * cosine_t_mat(1,4)     &
                                      +                    &
                       imag_v(i,2) * cosine_t_mat(2,4)     &
                                      +                    &
                       imag_v(i,3) * cosine_t_mat(3,4)     &
                                      +                    &
                       imag_v(i,4) * cosine_t_mat(4,4)     &
                                      +                    &
                       imag_v(i,5) * cosine_t_mat(5,4)     &
                                      +                    &
                       imag_v(i,6) * cosine_t_mat(6,4)     &
                                      +                    &
                       imag_v(i,7) * cosine_t_mat(7,4)     &
                                      +                    &
                       imag_v(i,8) * cosine_t_mat(8,4)     &
                                      +                    &
                       real_v(i,1) * sine_t_mat(1,4)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,4)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,4)       &
                                      +                    &
                       real_v(i,4) * sine_t_mat(4,4)       &
                                      +                    &
                       real_v(i,5) * sine_t_mat(5,4)       &
                                      +                    &
                       real_v(i,6) * sine_t_mat(6,4)       &
                                      +                    &
                       real_v(i,7) * sine_t_mat(7,4)       &
                                      +                    &
                       real_v(i,8) * sine_t_mat(8,4) 
     imag_v_scr(i,5) = imag_v(i,1) * cosine_t_mat(1,5)     &
                                      +                    &
                       imag_v(i,2) * cosine_t_mat(2,5)     &
                                      +                    &
                       imag_v(i,3) * cosine_t_mat(3,5)     &
                                      +                    &
                       imag_v(i,4) * cosine_t_mat(4,5)     &
                                      +                    &
                       imag_v(i,5) * cosine_t_mat(5,5)     &
                                      +                    &
                       imag_v(i,6) * cosine_t_mat(6,5)     &
                                      +                    &
                       imag_v(i,7) * cosine_t_mat(7,5)     &
                                      +                    &
                       imag_v(i,8) * cosine_t_mat(8,5)     &
                                      +                    &
                       real_v(i,1) * sine_t_mat(1,5)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,5)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,5)       &
                                      +                    &
                       real_v(i,4) * sine_t_mat(4,5)       &
                                      +                    &
                       real_v(i,5) * sine_t_mat(5,5)       &
                                      +                    &
                       real_v(i,6) * sine_t_mat(6,5)       &
                                      +                    &
                       real_v(i,7) * sine_t_mat(7,5)       &
                                      +                    &
                       real_v(i,8) * sine_t_mat(8,5) 
     imag_v_scr(i,6) = imag_v(i,1) * cosine_t_mat(1,6)     &
                                      +                    &
                       imag_v(i,2) * cosine_t_mat(2,6)     &
                                      +                    &
                       imag_v(i,3) * cosine_t_mat(3,6)     &
                                      +                    &
                       imag_v(i,4) * cosine_t_mat(4,6)     &
                                      +                    &
                       imag_v(i,5) * cosine_t_mat(5,6)     &
                                      +                    &
                       imag_v(i,6) * cosine_t_mat(6,6)     &
                                      +                    &
                       imag_v(i,7) * cosine_t_mat(7,6)     &
                                      +                    &
                       imag_v(i,8) * cosine_t_mat(8,6)     &
                                      +                    &
                       real_v(i,1) * sine_t_mat(1,6)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,6)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,6)       &
                                      +                    &
                       real_v(i,4) * sine_t_mat(4,6)       &
                                      +                    &
                       real_v(i,5) * sine_t_mat(5,6)       &
                                      +                    &
                       real_v(i,6) * sine_t_mat(6,6)       &
                                      +                    &
                       real_v(i,7) * sine_t_mat(7,6)       &
                                      +                    &
                       real_v(i,8) * sine_t_mat(8,6) 
     imag_v_scr(i,7) = imag_v(i,1) * cosine_t_mat(1,7)     &
                                      +                    &
                       imag_v(i,2) * cosine_t_mat(2,7)     &
                                      +                    &
                       imag_v(i,3) * cosine_t_mat(3,7)     &
                                      +                    &
                       imag_v(i,4) * cosine_t_mat(4,7)     &
                                      +                    &
                       imag_v(i,5) * cosine_t_mat(5,7)     &
                                      +                    &
                       imag_v(i,6) * cosine_t_mat(6,7)     &
                                      +                    &
                       imag_v(i,7) * cosine_t_mat(7,7)     &
                                      +                    &
                       imag_v(i,8) * cosine_t_mat(8,7)     &
                                      +                    &
                       real_v(i,1) * sine_t_mat(1,7)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,7)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,7)       &
                                      +                    &
                       real_v(i,4) * sine_t_mat(4,7)       &
                                      +                    &
                       real_v(i,5) * sine_t_mat(5,7)       &
                                      +                    &
                       real_v(i,6) * sine_t_mat(6,7)       &
                                      +                    &
                       real_v(i,7) * sine_t_mat(7,7)       &
                                      +                    &
                       real_v(i,8) * sine_t_mat(8,7) 
     imag_v_scr(i,8) = imag_v(i,1) * cosine_t_mat(1,8)     &
                                      +                    &
                       imag_v(i,2) * cosine_t_mat(2,8)     &
                                      +                    &
                       imag_v(i,3) * cosine_t_mat(3,8)     &
                                      +                    &
                       imag_v(i,4) * cosine_t_mat(4,8)     &
                                      +                    &
                       imag_v(i,5) * cosine_t_mat(5,8)     &
                                      +                    &
                       imag_v(i,6) * cosine_t_mat(6,8)     &
                                      +                    &
                       imag_v(i,7) * cosine_t_mat(7,8)     &
                                      +                    &
                       imag_v(i,8) * cosine_t_mat(8,8)     &
                                      +                    &
                       real_v(i,1) * sine_t_mat(1,8)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,8)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,8)       &
                                      +                    &
                       real_v(i,4) * sine_t_mat(4,8)       &
                                      +                    &
                       real_v(i,5) * sine_t_mat(5,8)       &
                                      +                    &
                       real_v(i,6) * sine_t_mat(6,8)       &
                                      +                    &
                       real_v(i,7) * sine_t_mat(7,8)       &
                                      +                    &
                       real_v(i,8) * sine_t_mat(8,8) 
  END DO
!
!       Copy the temporary vector back to the input vector.
!
  real_v(:,1:8) = real_v_scr(:,1:8) 
  imag_v(:,1:8) = imag_v_scr(:,1:8) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_m_8
