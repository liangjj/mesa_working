!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE v_v_m
!**begin prologue     v_v_m
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       v_v_m
!

  CONTAINS
!*deck v_v_m_gen
!***begin prologue     v_v_m_gen     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multipiles
!***                   for the general FEDVR HAmiltonian. 
!***                   The routine is not needed in 1D problems. 
!***                   In 2D, the i index runs over the number
!                      of points in the y coordinate while in 3D, it runs over
!                      the product of the number of points in the z and y coordinates.
!
!***references
!***routines called    ebcxx, ambcxx, apbcxx
!***end prologue   
!
  SUBROUTINE v_v_m_gen(real_v,imag_v,           &
                       real_v_scr,imag_v_scr,   &
                       cosine_t_mat,sine_t_mat, & 
                       ni,nj,nk)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, nk
  INTEGER                                  :: i, j, k
  REAL*8, DIMENSION(ni,nj)                 :: real_v, imag_v
  REAL*8, DIMENSION(ni,nj)                 :: real_v_scr, imag_v_scr
  REAL*8, DIMENSION(nk,nk)                 :: cosine_t_mat, sine_t_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
!
!                 Real(V_scr) = Real(V) * cos_t_mat 
!

!  CALL ebcxx(real_v_scr,real_v,cosine_t_mat,ni,nk,nk,ni,ni,nk)
!
!                 Real(V_scr) = Real(V_scr) - Imag(V) * sin_t_mat
!
!  CALL ambcxx(real_v_scr,imag_v,sine_t_mat,ni,nk,nk,ni,ni,nk)
!
!                 Imag(V_scr) = Imag(V) * cos_t_mat
!
!  CALL ebcxx(imag_v_scr,imag_v,cosine_t_mat,ni,nk,nk,ni,ni,nk)            
!
!                 Imag(V_scr) = Imag(V_scr) + Real(V) * sine_t_mat 
!
!  CALL apbcxx(imag_v_scr,real_v,sine_t_mat,ni,nk,nk,ni,ni,nk)      
!
!       Copy the temporary vector back to the input vector.
!
   real_v_scr(1:ni,1:nk) = 0.d0
   imag_v_scr(1:ni,1:nk) = 0.d0
   DO k=1,nk
      DO i=1,ni
         DO j=1,nk
            real_v_scr(i,j) = real_v_scr(i,j) + cosine_t_mat(k,j) * real_v(i,k) &
                                              -   sine_t_mat(k,j) * imag_v(i,k)
            imag_v_scr(i,j) = imag_v_scr(i,j) + cosine_t_mat(k,j) * imag_v(i,k) &
                                              +   sine_t_mat(k,j) * real_v(i,k)
         END DO
      END DO
   END DO
  real_v(:,1:nk) = real_v_scr(:,1:nk) 
  imag_v(:,1:nk) = imag_v_scr(:,1:nk) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE v_v_m_gen
!
!
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
!
!
!deck v_v_m_3
!***begin prologue     v_v_m_3     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_v_m_3
!
  SUBROUTINE v_v_m_3(real_v,imag_v,real_v_scr,imag_v_scr, &
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
  DO i=1,ni
     real_v_scr(i,1) = real_v(i,1) * cosine_t_mat(1,1)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,1)     &
                                      +                    &
                       real_v(i,3) * cosine_t_mat(3,1)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,1)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,1)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,1)    
     real_v_scr(i,2) = real_v(i,1) * cosine_t_mat(1,2)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,2)     &
                                      +                    &
                       real_v(i,3) * cosine_t_mat(3,2)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,2)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,2)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,2)    
     real_v_scr(i,3) = real_v(i,1) * cosine_t_mat(1,3)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,3)     &
                                      +                    &
                       real_v(i,3) * cosine_t_mat(3,3)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,3)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,3)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,3)    

     imag_v_scr(i,1) = imag_v(i,1) * cosine_t_mat(1,1)     &
                                      +                    &
                       imag_v(i,2) * cosine_t_mat(2,1)     &
                                      +                    &
                       imag_v(i,3) * cosine_t_mat(3,1)     &
                                      +                    &
                       real_v(i,1) * sine_t_mat(1,1)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,1)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,1) 
     imag_v_scr(i,2) = imag_v(i,1) * cosine_t_mat(1,2)     &
                                      +                    &
                       imag_v(i,2) * cosine_t_mat(2,2)     &
                                      +                    &
                       imag_v(i,3) * cosine_t_mat(3,2)     &
                                      +                    &
                       real_v(i,1) * sine_t_mat(1,2)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,2)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,2) 
     imag_v_scr(i,3) = imag_v(i,1) * cosine_t_mat(1,3)     &
                                      +                    &
                       imag_v(i,2) * cosine_t_mat(2,3)     &
                                      +                    &
                       imag_v(i,3) * cosine_t_mat(3,3)     &
                                      +                    &
                       real_v(i,1) * sine_t_mat(1,3)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,3)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,3) 

  END DO
!
!       Copy the temporary vector back to the input vector.
!
  real_v(:,1:3) = real_v_scr(:,1:3) 
  imag_v(:,1:3) = imag_v_scr(:,1:3) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_m_3
!
!
!deck v_v_m_4
!***begin prologue     v_v_m_4     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Computes the vector propagator matrix multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_v_m_4
!
  SUBROUTINE v_v_m_4(real_v,imag_v,real_v_scr,imag_v_scr, &
                     cosine_t_mat,sine_t_mat,ni,nj)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: real_v, imag_v
  REAL*8, DIMENSION(ni,nj)                 :: real_v_scr, imag_v_scr
  REAL*8, DIMENSION(4,4)                   :: cosine_t_mat
  REAL*8, DIMENSION(4,4)                   :: sine_t_mat
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
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,1)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,1)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,1)       &
                                      -                    &
                       imag_v(i,4) * sine_t_mat(4,1)    
     real_v_scr(i,2) = real_v(i,1) * cosine_t_mat(1,2)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,2)     &
                                      +                    &
                       real_v(i,3) * cosine_t_mat(3,2)     &
                                      +                    &
                       real_v(i,4) * cosine_t_mat(4,2)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,2)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,2)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,2)       &
                                      -                    & 
                       imag_v(i,4) * sine_t_mat(4,2)    
     real_v_scr(i,3) = real_v(i,1) * cosine_t_mat(1,3)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,3)     &
                                      +                    &
                       real_v(i,3) * cosine_t_mat(3,3)     &
                                      +                    &
                       real_v(i,4) * cosine_t_mat(4,3)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,3)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,3)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,3)       &
                                      -                    &
                       imag_v(i,4) * sine_t_mat(4,3)    
     real_v_scr(i,4) = real_v(i,1) * cosine_t_mat(1,4)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,4)     &
                                      +                    &
                       real_v(i,3) * cosine_t_mat(3,4)     &
                                      +                    &
                       real_v(i,4) * cosine_t_mat(4,4)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,4)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,4)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,4)       &
                                      -                    &
                       imag_v(i,4) * sine_t_mat(4,4)    

     imag_v_scr(i,1) = imag_v(i,1) * cosine_t_mat(1,1)     &
                                      +                    &
                       imag_v(i,2) * cosine_t_mat(2,1)     &
                                      +                    &
                       imag_v(i,3) * cosine_t_mat(3,1)     &
                                      +                    &
                       imag_v(i,4) * cosine_t_mat(4,1)     &
                                      +                    &
                       real_v(i,1) * sine_t_mat(1,1)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,1)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,1)       &
                                      +                    &
                       real_v(i,4) * sine_t_mat(4,1) 
     imag_v_scr(i,2) = imag_v(i,1) * cosine_t_mat(1,2)     &
                                      +                    &
                       imag_v(i,2) * cosine_t_mat(2,2)     &
                                      +                    &
                       imag_v(i,3) * cosine_t_mat(3,2)     &
                                      +                    &
                       imag_v(i,4) * cosine_t_mat(4,2)     &
                                      +                    &
                       real_v(i,1) * sine_t_mat(1,2)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,2)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,2)       &
                                      +                    &
                       real_v(i,4) * sine_t_mat(4,2) 
     imag_v_scr(i,3) = imag_v(i,1) * cosine_t_mat(1,3)     &
                                      +                    &
                       imag_v(i,2) * cosine_t_mat(2,3)     &
                                      +                    &
                       imag_v(i,3) * cosine_t_mat(3,3)     &
                                      +                    &
                       imag_v(i,4) * cosine_t_mat(4,3)     &
                                      +                    &
                       real_v(i,1) * sine_t_mat(1,3)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,3)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,3)       &
                                      +                    &
                       real_v(i,4) * sine_t_mat(4,3)       
     imag_v_scr(i,4) = imag_v(i,1) * cosine_t_mat(1,4)     &
                                      +                    &
                       imag_v(i,2) * cosine_t_mat(2,4)     &
                                      +                    &
                       imag_v(i,3) * cosine_t_mat(3,4)     &
                                      +                    &
                       imag_v(i,4) * cosine_t_mat(4,4)     &
                                      +                    &
                       real_v(i,1) * sine_t_mat(1,4)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,4)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,4)       &
                                      +                    &
                       real_v(i,4) * sine_t_mat(4,4)       

  END DO
!
!       Copy the temporary vector back to the input vector.
!
  real_v(:,1:4) = real_v_scr(:,1:4) 
  imag_v(:,1:4) = imag_v_scr(:,1:4) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_m_4
!
!
!deck v_v_m_5
!***begin prologue     v_v_m_5     
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
!***end prologue       v_v_m_5
!
  SUBROUTINE v_v_m_5(real_v,imag_v,real_v_scr,imag_v_scr, &
                     cosine_t_mat,sine_t_mat,ni,nj)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: real_v, imag_v
  REAL*8, DIMENSION(ni,nj)                 :: real_v_scr, imag_v_scr
  REAL*8, DIMENSION(5,5)                   :: cosine_t_mat
  REAL*8, DIMENSION(5,5)                   :: sine_t_mat
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
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,1)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,1)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,1)       &
                                      -                    &
                       imag_v(i,4) * sine_t_mat(4,1)       &
                                      -                    &
                       imag_v(i,5) * sine_t_mat(5,1)    
     real_v_scr(i,2) = real_v(i,1) * cosine_t_mat(1,2)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,2)     &
                                      +                    &
                       real_v(i,3) * cosine_t_mat(3,2)     &
                                      +                    &
                       real_v(i,4) * cosine_t_mat(4,2)     &
                                      +                    &
                       real_v(i,5) * cosine_t_mat(5,2)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,2)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,2)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,2)       &
                                      -                    & 
                       imag_v(i,4) * sine_t_mat(4,2)       &
                                      -                    & 
                       imag_v(i,5) * sine_t_mat(5,2)    
     real_v_scr(i,3) = real_v(i,1) * cosine_t_mat(1,3)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,3)     &
                                      +                    &
                       real_v(i,3) * cosine_t_mat(3,3)     &
                                      +                    &
                       real_v(i,4) * cosine_t_mat(4,3)     &
                                      +                    &
                       real_v(i,5) * cosine_t_mat(5,3)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,3)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,3)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,3)       &
                                      -                    &
                       imag_v(i,4) * sine_t_mat(4,3)       &
                                      -                    &
                       imag_v(i,5) * sine_t_mat(5,3)    
     real_v_scr(i,4) = real_v(i,1) * cosine_t_mat(1,4)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,4)     &
                                      +                    &
                       real_v(i,3) * cosine_t_mat(3,4)     &
                                      +                    &
                       real_v(i,4) * cosine_t_mat(4,4)     &
                                      +                    &
                       real_v(i,5) * cosine_t_mat(5,4)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,4)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,4)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,4)       &
                                      -                    &
                       imag_v(i,4) * sine_t_mat(4,4)       &
                                      -                    &
                       imag_v(i,5) * sine_t_mat(5,4)    
     real_v_scr(i,5) = real_v(i,1) * cosine_t_mat(1,5)     &
                                      +                    &
                       real_v(i,2) * cosine_t_mat(2,5)     &
                                      +                    &
                       real_v(i,3) * cosine_t_mat(3,5)     &
                                      +                    &
                       real_v(i,4) * cosine_t_mat(4,5)     &
                                      +                    &
                       real_v(i,5) * cosine_t_mat(5,5)     &
                                      -                    &
                       imag_v(i,1) * sine_t_mat(1,5)       &
                                      -                    &
                       imag_v(i,2) * sine_t_mat(2,5)       &
                                      -                    &
                       imag_v(i,3) * sine_t_mat(3,5)       &
                                      -                    &
                       imag_v(i,4) * sine_t_mat(4,5)       &
                                      -                    &
                       imag_v(i,5) * sine_t_mat(5,5)    

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
                       real_v(i,1) * sine_t_mat(1,1)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,1)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,1)       &
                                      +                    &
                       real_v(i,4) * sine_t_mat(4,1)       &
                                      +                    &
                       real_v(i,5) * sine_t_mat(5,1) 
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
                       real_v(i,1) * sine_t_mat(1,2)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,2)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,2)       &
                                      +                    &
                       real_v(i,4) * sine_t_mat(4,2)       &
                                      +                    &
                       real_v(i,5) * sine_t_mat(5,2) 
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
                       real_v(i,1) * sine_t_mat(1,3)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,3)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,3)       &
                                      +                    &
                       real_v(i,4) * sine_t_mat(4,3)       &
                                      +                    &
                       real_v(i,5) * sine_t_mat(5,3)       
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
                       real_v(i,1) * sine_t_mat(1,4)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,4)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,4)       &
                                      +                    &
                       real_v(i,4) * sine_t_mat(4,4)       &
                                      +                    &
                       real_v(i,5) * sine_t_mat(5,4)       
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
                       real_v(i,1) * sine_t_mat(1,5)       &
                                      +                    &
                       real_v(i,2) * sine_t_mat(2,5)       &
                                      +                    &
                       real_v(i,3) * sine_t_mat(3,5)       &
                                      +                    &
                       real_v(i,4) * sine_t_mat(4,5)       &
                                      +                    &
                       real_v(i,5) * sine_t_mat(5,5)       

  END DO
!
!       Copy the temporary vector back to the input vector.
!
  real_v(:,1:5) = real_v_scr(:,1:5) 
  imag_v(:,1:5) = imag_v_scr(:,1:5) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_m_5
!
!
!deck v_v_m_6
!***begin prologue     v_v_m_6     
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
!***end prologue       v_v_m_6
!
  SUBROUTINE v_v_m_6(real_v,imag_v,real_v_scr,imag_v_scr, &
                     cosine_t_mat,sine_t_mat,ni,nj)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: real_v, imag_v
  REAL*8, DIMENSION(ni,nj)                 :: real_v_scr, imag_v_scr
  REAL*8, DIMENSION(6,6)                   :: cosine_t_mat
  REAL*8, DIMENSION(6,6)                   :: sine_t_mat
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
                       imag_v(i,6) * sine_t_mat(6,1)    
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
                       imag_v(i,6) * sine_t_mat(6,2)    
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
                       imag_v(i,6) * sine_t_mat(6,3)    
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
                       imag_v(i,6) * sine_t_mat(6,4)    
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
                       imag_v(i,6) * sine_t_mat(6,5)    
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
                       imag_v(i,6) * sine_t_mat(6,6)    

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
                       real_v(i,6) * sine_t_mat(6,1) 
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
                       real_v(i,6) * sine_t_mat(6,2) 
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
                       real_v(i,6) * sine_t_mat(6,3) 
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
                       real_v(i,6) * sine_t_mat(6,4) 
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
                       real_v(i,6) * sine_t_mat(6,5) 
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
                       real_v(i,6) * sine_t_mat(6,6) 
  END DO
!
!       Copy the temporary vector back to the input vector.
!
  real_v(:,1:6) = real_v_scr(:,1:6) 
  imag_v(:,1:6) = imag_v_scr(:,1:6) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_m_6
!
!
!deck v_v_m_7
!***begin prologue     v_v_m_7     
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
!***end prologue       v_v_m_7
!
  SUBROUTINE v_v_m_7(real_v,imag_v,real_v_scr,imag_v_scr, &
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
                       imag_v(i,7) * sine_t_mat(7,1)       
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
                       imag_v(i,7) * sine_t_mat(7,2)       
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
                       imag_v(i,7) * sine_t_mat(7,3)       
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
                       imag_v(i,7) * sine_t_mat(7,4)       
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
                       imag_v(i,7) * sine_t_mat(7,5)       
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
                       imag_v(i,7) * sine_t_mat(7,6)       
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
                       imag_v(i,7) * sine_t_mat(7,7)       
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
                       real_v(i,7) * sine_t_mat(7,1)       
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
                       real_v(i,7) * sine_t_mat(7,2)       
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
                       real_v(i,7) * sine_t_mat(7,3)       
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
                       real_v(i,7) * sine_t_mat(7,4)       
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
                       real_v(i,7) * sine_t_mat(7,5)       
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
                       real_v(i,7) * sine_t_mat(7,6)       
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
                       real_v(i,7) * sine_t_mat(7,7)       
  END DO
!
!       Copy the temporary vector back to the input vector.
!
  real_v(:,1:7) = real_v_scr(:,1:7) 
  imag_v(:,1:7) = imag_v_scr(:,1:7) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_m_7
!
!
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
!
!
!deck v_v_m_9
!***begin prologue     v_v_m_9     
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
!***end prologue       v_v_m_9
!
  SUBROUTINE v_v_m_9(real_v,imag_v,real_v_scr,imag_v_scr, &
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
  DO i=1,ni
     real_v_scr(i,1) = real_v(i,1)  * cosine_t_mat(1,1)    &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,1)    &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,1)    &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,1)    &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,1)    &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,1)    &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,1)    &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,1)    &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,1)    &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,1)      &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,1)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,1)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,1)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,1)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,1)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,1)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,1)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,1)      
     real_v_scr(i,2) = real_v(i,1)  * cosine_t_mat(1,2)    &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,2)    &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,2)    &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,2)    &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,2)    &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,2)    &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,2)    &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,2)    &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,2)    &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,2)      &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,2)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,2)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,2)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,2)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,2)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,2)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,2)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,2)      

     real_v_scr(i,3) = real_v(i,1)  * cosine_t_mat(1,3)    &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,3)    &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,3)    &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,3)    &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,3)    &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,3)    &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,3)    &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,3)    &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,3)    &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,3)      &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,3)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,3)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,3)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,3)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,3)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,3)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,3)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,3)      
     real_v_scr(i,4) = real_v(i,1)  * cosine_t_mat(1,4)    &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,4)    &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,4)    &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,4)    &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,4)    &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,4)    &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,4)    &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,4)    &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,4)    &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,4)      &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,4)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,4)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,4)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,4)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,4)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,4)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,4)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,4)      
     real_v_scr(i,5) = real_v(i,1)  * cosine_t_mat(1,5)    &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,5)    &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,5)    &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,5)    &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,5)    &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,5)    &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,5)    &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,5)    &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,5)    &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,5)      &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,5)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,5)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,5)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,5)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,5)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,5)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,5)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,5)      
     real_v_scr(i,6) = real_v(i,1)  * cosine_t_mat(1,6)    &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,6)    &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,6)    &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,6)    &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,6)    &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,6)    &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,6)    &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,6)    &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,6)    &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,6)      &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,6)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,6)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,6)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,6)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,6)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,6)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,6)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,6)      
     real_v_scr(i,7) = real_v(i,1)  * cosine_t_mat(1,7)    &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,7)    &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,7)    &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,7)    &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,7)    &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,7)    &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,7)    &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,7)    &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,7)    &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,7)      &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,7)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,7)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,7)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,7)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,7)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,7)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,7)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,7)      
     real_v_scr(i,8) = real_v(i,1)  * cosine_t_mat(1,8)    &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,8)    &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,8)    &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,8)    &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,8)    &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,8)    &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,8)    &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,8)    &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,8)    &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,8)      &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,8)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,8)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,8)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,8)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,8)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,8)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,8)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,8)      
     real_v_scr(i,9) = real_v(i,1)  * cosine_t_mat(1,9)    &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,9)    &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,9)    &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,9)    &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,9)    &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,9)    &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,9)    &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,9)    &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,9)    &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,9)      &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,9)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,9)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,9)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,9)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,9)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,9)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,9)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,9)      
     imag_v_scr(i,1) = imag_v(i,1)  * cosine_t_mat(1,1)    &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,1)    &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,1)    &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,1)    &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,1)    &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,1)    &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,1)    &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,1)    &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,1)    &
                                    +                      &
                       real_v(i,1)  * sine_t_mat(1,1)      &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,1)      &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,1)      &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,1)      &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,1)      &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,1)      &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,1)      &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,1)      &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,1)      
     imag_v_scr(i,2) = imag_v(i,1)  * cosine_t_mat(1,2)    &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,2)    &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,2)    &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,2)    &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,2)    &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,2)    &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,2)    &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,2)    &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,2)    &
                                    +                      &
                       real_v(i,1)  * sine_t_mat(1,2)      &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,2)      &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,2)      &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,2)      &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,2)      &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,2)      &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,2)      &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,2)      &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,2)      
     imag_v_scr(i,3) = imag_v(i,1)  * cosine_t_mat(1,3)    &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,3)    &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,3)    &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,3)    &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,3)    &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,3)    &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,3)    &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,3)    &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,3)    &
                                    +                      &
                       real_v(i,1)  * sine_t_mat(1,3)      &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,3)      &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,3)      &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,3)      &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,3)      &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,3)      &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,3)      &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,3)      &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,3)      
     imag_v_scr(i,4) = imag_v(i,1)  * cosine_t_mat(1,4)    &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,4)    &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,4)    &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,4)    &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,4)    &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,4)    &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,4)    &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,4)    &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,4)    &
                                    +                      &
                       real_v(i,1)  * sine_t_mat(1,4)      &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,4)      &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,4)      &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,4)      &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,4)      &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,4)      &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,4)      &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,4)      &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,4)      
     imag_v_scr(i,5) = imag_v(i,1)  * cosine_t_mat(1,5)    &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,5)    &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,5)    &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,5)    &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,5)    &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,5)    &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,5)    &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,5)    &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,5)    &
                                    +                      &
                       real_v(i,1)  * sine_t_mat(1,5)      &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,5)      &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,5)      &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,5)      &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,5)      &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,5)      &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,5)      &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,5)      &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,5)      
     imag_v_scr(i,6) = imag_v(i,1)  * cosine_t_mat(1,6)    &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,6)    &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,6)    &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,6)    &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,6)    &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,6)    &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,6)    &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,6)    &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,6)    &
                                    +                      &
                       real_v(i,1)  * sine_t_mat(1,6)      &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,6)      &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,6)      &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,6)      &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,6)      &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,6)      &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,6)      &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,6)      &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,6)      
     imag_v_scr(i,7) = imag_v(i,1)  * cosine_t_mat(1,7)    &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,7)    &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,7)    &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,7)    &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,7)    &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,7)    &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,7)    &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,7)    &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,7)    &
                                    +                      &
                       real_v(i,1)  * sine_t_mat(1,7)      &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,7)      &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,7)      &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,7)      &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,7)      &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,7)      &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,7)      &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,7)      &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,7)      
     imag_v_scr(i,8) = imag_v(i,1)  * cosine_t_mat(1,8)    &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,8)    &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,8)    &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,8)    &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,8)    &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,8)    &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,8)    &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,8)    &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,8)    &
                                    +                      &
                       real_v(i,1)  * sine_t_mat(1,8)      &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,8)      &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,8)      &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,8)      &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,8)      &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,8)      &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,8)      &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,8)      &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,8)      
     imag_v_scr(i,9) = imag_v(i,1)  * cosine_t_mat(1,9)    &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,9)    &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,9)    &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,9)    &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,9)    &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,9)    &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,9)    &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,9)    &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,9)    &
                                    +                      &
                       real_v(i,1)  * sine_t_mat(1,9)      &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,9)      &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,9)      &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,9)      &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,9)      &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,9)      &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,9)      &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,9)      &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,9)      
  END DO
!
!       Copy the temporary vector back to the input vector.
!
  real_v(:,1:9) = real_v_scr(:,1:9) 
  imag_v(:,1:9) = imag_v_scr(:,1:9) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_m_9
!
!
!deck v_v_m_10
!***begin prologue     v_v_m_10     
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
!***end prologue       v_v_m_10
!
  SUBROUTINE v_v_m_10(real_v,imag_v,real_v_scr,imag_v_scr, &
                     cosine_t_mat,sine_t_mat,ni,nj)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: real_v, imag_v
  REAL*8, DIMENSION(ni,nj)                 :: real_v_scr, imag_v_scr
  REAL*8, DIMENSION(10,10)                 :: cosine_t_mat
  REAL*8, DIMENSION(10,10)                 :: sine_t_mat
  INTEGER                                  :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  DO i=1,ni
     real_v_scr(i,1) = real_v(i,1)  * cosine_t_mat(1,1)    &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,1)    &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,1)    &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,1)    &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,1)    &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,1)    &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,1)    &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,1)    &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,1)    &
                                    +                      &
                       real_v(i,10) * cosine_t_mat(10,1)   &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,1)      &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,1)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,1)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,1)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,1)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,1)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,1)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,1)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,1)      &
                                    -                      &
                       imag_v(i,10) * sine_t_mat(10,1)    
     real_v_scr(i,2) = real_v(i,1)  * cosine_t_mat(1,2)    &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,2)    &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,2)    &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,2)    &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,2)    &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,2)    &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,2)    &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,2)    &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,2)    &
                                    +                      &
                       real_v(i,10) * cosine_t_mat(10,2)   &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,2)      &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,2)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,2)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,2)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,2)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,2)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,2)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,2)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,2)      &
                                    -                      &
                       imag_v(i,10) * sine_t_mat(10,2)     &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,2)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,2)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,2)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,2)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,2)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,2)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,2)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,2)      &
                                    -                      &
                       imag_v(i,10) * sine_t_mat(10,2)    
     real_v_scr(i,3) = real_v(i,1)  * cosine_t_mat(1,3)    &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,3)    &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,3)    &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,3)    &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,3)    &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,3)    &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,3)    &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,3)    &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,3)    &
                                    +                      &
                       real_v(i,10) * cosine_t_mat(10,3)   &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,3)      &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,3)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,3)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,3)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,3)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,3)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,3)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,3)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,3)      &
                                    -                      &
                       imag_v(i,10) * sine_t_mat(10,3)    
     real_v_scr(i,4) = real_v(i,1)  * cosine_t_mat(1,4)    &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,4)    &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,4)    &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,4)    &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,4)    &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,4)    &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,4)    &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,4)    &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,4)    &
                                    +                      &
                       real_v(i,10) * cosine_t_mat(10,4)   &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,4)      &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,4)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,4)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,4)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,4)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,4)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,4)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,4)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,4)      &
                                    -                      &
                       imag_v(i,10) * sine_t_mat(10,4)    
     real_v_scr(i,5) = real_v(i,1)  * cosine_t_mat(1,5)    &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,5)    &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,5)    &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,5)    &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,5)    &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,5)    &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,5)    &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,5)    &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,5)    &
                                    +                      &
                       real_v(i,10) * cosine_t_mat(10,5)   &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,5)      &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,5)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,5)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,5)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,5)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,5)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,5)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,5)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,5)      &
                                    -                      &
                       imag_v(i,10) * sine_t_mat(10,5)    
     real_v_scr(i,6) = real_v(i,1)  * cosine_t_mat(1,6)    &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,6)    &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,6)    &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,6)    &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,6)    &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,6)    &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,6)    &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,6)    &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,6)    &
                                    +                      &
                       real_v(i,10) * cosine_t_mat(10,6)   &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,6)      &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,6)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,6)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,6)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,6)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,6)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,6)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,6)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,6)      &
                                    -                      &
                       imag_v(i,10) * sine_t_mat(10,6)    
     real_v_scr(i,7) = real_v(i,1)  * cosine_t_mat(1,7)    &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,7)    &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,7)    &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,7)    &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,7)    &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,7)    &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,7)    &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,7)    &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,7)    &
                                    +                      &
                       real_v(i,10) * cosine_t_mat(10,7)   &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,7)      &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,7)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,7)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,7)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,7)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,7)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,7)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,7)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,7)      &
                                    -                      &
                       imag_v(i,10) * sine_t_mat(10,7)    
     real_v_scr(i,8) = real_v(i,1)  * cosine_t_mat(1,8)    &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,8)    &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,8)    &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,8)    &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,8)    &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,8)    &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,8)    &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,8)    &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,8)    &
                                    +                      &
                       real_v(i,10) * cosine_t_mat(10,8)   &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,8)      &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,8)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,8)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,8)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,8)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,8)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,8)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,8)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,8)      &
                                    -                      &
                       imag_v(i,10) * sine_t_mat(10,8)    
     real_v_scr(i,9) = real_v(i,1)  * cosine_t_mat(1,9)    &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,9)    &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,9)    &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,9)    &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,9)    &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,9)    &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,9)    &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,9)    &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,9)    &
                                    +                      &
                       real_v(i,10) * cosine_t_mat(10,9)   &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,9)      &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,9)      &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,9)      &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,9)      &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,9)      &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,9)      &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,9)      &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,9)      &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,9)      &
                                    -                      &
                       imag_v(i,10) * sine_t_mat(10,9)    
     real_v_scr(i,10) = real_v(i,1) * cosine_t_mat(1,10)   &
                                    +                      &
                       real_v(i,2)  * cosine_t_mat(2,10)   &
                                    +                      &
                       real_v(i,3)  * cosine_t_mat(3,10)   &
                                    +                      &
                       real_v(i,4)  * cosine_t_mat(4,10)   &
                                    +                      &
                       real_v(i,5)  * cosine_t_mat(5,10)   &
                                    +                      &
                       real_v(i,6)  * cosine_t_mat(6,10)   &
                                    +                      &
                       real_v(i,7)  * cosine_t_mat(7,10)   &
                                    +                      &
                       real_v(i,8)  * cosine_t_mat(8,10)   &
                                    +                      &
                       real_v(i,9)  * cosine_t_mat(9,10)   &
                                    +                      &
                       real_v(i,10) * cosine_t_mat(10,10)  &
                                    -                      &
                       imag_v(i,1)  * sine_t_mat(1,10)     &
                                    -                      &
                       imag_v(i,2)  * sine_t_mat(2,10)     &
                                    -                      &
                       imag_v(i,3)  * sine_t_mat(3,10)     &
                                    -                      &
                       imag_v(i,4)  * sine_t_mat(4,10)     &
                                    -                      &
                       imag_v(i,5)  * sine_t_mat(5,10)     &
                                    -                      &
                       imag_v(i,6)  * sine_t_mat(6,10)     &
                                    -                      &
                       imag_v(i,7)  * sine_t_mat(7,10)     &
                                    -                      &
                       imag_v(i,8)  * sine_t_mat(8,10)     &
                                    -                      &
                       imag_v(i,9)  * sine_t_mat(9,10)     &
                                    -                      &
                       imag_v(i,10) * sine_t_mat(10,10)    
     imag_v_scr(i,1) = imag_v(i,1)  * cosine_t_mat(1,1)    &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,1)    &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,1)    &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,1)    &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,1)    &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,1)    &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,1)    &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,1)    &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,1)    &
                                    +                      &
                       imag_v(i,10) * cosine_t_mat(10,1)   &
                                    +                      &
                       real_v(i,1)  * sine_t_mat(1,1)      &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,1)      &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,1)      &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,1)      &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,1)      &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,1)      &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,1)      &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,1)      &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,1)      &
                                    +                      &
                       real_v(i,10) * sine_t_mat(10,1) 
     imag_v_scr(i,2) = imag_v(i,1)  * cosine_t_mat(1,2)    &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,2)    &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,2)    &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,2)    &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,2)    &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,2)    &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,2)    &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,2)    &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,2)    &
                                    +                      &
                       imag_v(i,10) * cosine_t_mat(10,2)   &
                                    +                      &
                       real_v(i,1)  * sine_t_mat(1,2)      &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,2)      &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,2)      &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,2)      &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,2)      &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,2)      &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,2)      &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,2)      &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,2)      &
                                    +                      &
                       real_v(i,10) * sine_t_mat(10,2) 
     imag_v_scr(i,3) = imag_v(i,1)  * cosine_t_mat(1,3)    &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,3)    &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,3)    &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,3)    &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,3)    &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,3)    &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,3)    &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,3)    &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,3)    &
                                    +                      &
                       imag_v(i,10) * cosine_t_mat(10,3)   &
                                    +                      &
                       real_v(i,1)  * sine_t_mat(1,3)      &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,3)      &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,3)      &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,3)      &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,3)      &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,3)      &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,3)      &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,3)      &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,3)      &
                                    +                      &
                       real_v(i,10) * sine_t_mat(10,3) 
     imag_v_scr(i,4) = imag_v(i,1)  * cosine_t_mat(1,4)    &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,4)    &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,4)    &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,4)    &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,4)    &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,4)    &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,4)    &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,4)    &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,4)    &
                                    +                      &
                       imag_v(i,10) * cosine_t_mat(10,4)   &
                                    +                      &
                       real_v(i,1)  * sine_t_mat(1,4)      &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,4)      &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,4)      &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,4)      &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,4)      &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,4)      &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,4)      &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,4)      &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,4)      &
                                    +                      &
                       real_v(i,10) * sine_t_mat(10,4) 
     imag_v_scr(i,5) = imag_v(i,1)  * cosine_t_mat(1,5)    &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,5)    &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,5)    &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,5)    &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,5)    &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,5)    &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,5)    &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,5)    &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,5)    &
                                    +                      &
                       imag_v(i,10) * cosine_t_mat(10,5)   &
                                    +                      &
                       real_v(i,1)  * sine_t_mat(1,5)      &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,5)      &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,5)      &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,5)      &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,5)      &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,5)      &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,5)      &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,5)      &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,5)      &
                                    +                      &
                       real_v(i,10) * sine_t_mat(10,5) 
     imag_v_scr(i,6) = imag_v(i,1)  * cosine_t_mat(1,6)    &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,6)    &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,6)    &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,6)    &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,6)    &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,6)    &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,6)    &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,6)    &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,6)    &
                                    +                      &
                       imag_v(i,10) * cosine_t_mat(10,6)   &
                                    +                      &
                       real_v(i,1)  * sine_t_mat(1,6)      &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,6)      &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,6)      &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,6)      &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,6)      &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,6)      &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,6)      &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,6)      &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,6)      &
                                    +                      &
                       real_v(i,10) * sine_t_mat(10,6) 
     imag_v_scr(i,7) = imag_v(i,1)  * cosine_t_mat(1,7)    &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,7)    &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,7)    &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,7)    &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,7)    &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,7)    &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,7)    &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,7)    &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,7)    &
                                    +                      &
                       imag_v(i,10) * cosine_t_mat(10,7)   &
                                    +                      &
                       real_v(i,1)  * sine_t_mat(1,7)      &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,7)      &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,7)      &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,7)      &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,7)      &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,7)      &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,7)      &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,7)      &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,7)      &
                                    +                      &
                       real_v(i,10) * sine_t_mat(10,7) 
     imag_v_scr(i,8) = imag_v(i,1)  * cosine_t_mat(1,8)    &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,8)    &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,8)    &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,8)    &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,8)    &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,8)    &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,8)    &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,8)    &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,8)    &
                                    +                      &
                       imag_v(i,10) * cosine_t_mat(10,8)   &
                                    +                      &
                       real_v(i,1)  * sine_t_mat(1,8)      &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,8)      &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,8)      &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,8)      &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,8)      &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,8)      &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,8)      &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,8)      &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,8)      &
                                    +                      &
                       real_v(i,10) * sine_t_mat(10,8) 
     imag_v_scr(i,9) = imag_v(i,1)  * cosine_t_mat(1,9)    &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,9)    &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,9)    &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,9)    &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,9)    &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,9)    &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,9)    &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,9)    &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,9)    &
                                    +                      &
                       imag_v(i,10) * cosine_t_mat(10,9)   &
                                    +                      &
                       real_v(i,1)  * sine_t_mat(1,9)      &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,9)      &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,9)      &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,9)      &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,9)      &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,9)      &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,9)      &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,9)      &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,9)      &
                                    +                      &
                       real_v(i,10) * sine_t_mat(10,9) 
     imag_v_scr(i,10) = imag_v(i,1) * cosine_t_mat(1,10)   &
                                    +                      &
                       imag_v(i,2)  * cosine_t_mat(2,10)   &
                                    +                      &
                       imag_v(i,3)  * cosine_t_mat(3,10)   &
                                    +                      &
                       imag_v(i,4)  * cosine_t_mat(4,10)   &
                                    +                      &
                       imag_v(i,5)  * cosine_t_mat(5,10)   &
                                    +                      &
                       imag_v(i,6)  * cosine_t_mat(6,10)   &
                                    +                      &
                       imag_v(i,7)  * cosine_t_mat(7,10)   &
                                    +                      &
                       imag_v(i,8)  * cosine_t_mat(8,10)   &
                                    +                      &
                       imag_v(i,9)  * cosine_t_mat(9,10)   &
                                    +                      &
                       imag_v(i,10) * cosine_t_mat(10,10)  &
                                    +                      &
                       real_v(i,10) * sine_t_mat(1,10)     &
                                    +                      &
                       real_v(i,2)  * sine_t_mat(2,10)     &
                                    +                      &
                       real_v(i,3)  * sine_t_mat(3,10)     &
                                    +                      &
                       real_v(i,4)  * sine_t_mat(4,10)     &
                                    +                      &
                       real_v(i,5)  * sine_t_mat(5,10)     &
                                    +                      &
                       real_v(i,6)  * sine_t_mat(6,10)     &
                                    +                      &
                       real_v(i,7)  * sine_t_mat(7,10)     &
                                    +                      &
                       real_v(i,8)  * sine_t_mat(8,10)     &
                                    +                      &
                       real_v(i,9)  * sine_t_mat(9,10)     &
                                    +                      &
                       real_v(i,10) * sine_t_mat(10,10) 
  END DO
!
!       Copy the temporary vector back to the input vector.
!
  real_v(:,1:10) = real_v_scr(:,1:10) 
  imag_v(:,1:10) = imag_v_scr(:,1:10) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_m_10
!
END MODULE v_v_m
