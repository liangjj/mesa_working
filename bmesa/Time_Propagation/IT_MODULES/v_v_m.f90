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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  SUBROUTINE v_v_m_gen(v,          &
                       v_scr,      &
                       exp_t_mat,  & 
                       ni,nj,nk)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, nk
  INTEGER                                  :: i, j, k
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(nk,nk)                 :: exp_t_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
!
!                 Real(V_scr) = Real(V) * cos_t_mat 
!

!  CALL ebcxx(v_scr,v,cosine_t_mat,ni,nk,nk,ni,ni,nk)
!
!                 Real(V_scr) = Real(V_scr) - Imag(V) * sin_t_mat
!
!  CALL ambcxx(v_scr,sine_t_mat,ni,nk,nk,ni,ni,nk)
!
!                 Imag(V_scr) = Imag(V) * cos_t_mat
!
!  CALL ebcxx(imag_v_scr,cosine_t_mat,ni,nk,nk,ni,ni,nk)            
!
!                 Imag(V_scr) = Imag(V_scr) + Real(V) * sine_t_mat 
!
!  CALL apbcxx(imag_v_scr,v,sine_t_mat,ni,nk,nk,ni,ni,nk)      
!
!       Copy the temporary vector back to the input vector.
!
   v_scr(1:ni,1:nk) = 0.d0
   DO k=1,nk
      DO i=1,ni
         DO j=1,nk
            v_scr(i,j) = v_scr(i,j) + exp_t_mat(k,j) * v(i,k) 
         END DO
      END DO
   END DO
  v(:,1:nk) = v_scr(:,1:nk) 
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
  SUBROUTINE v_v_m_2(v,          &
                     v_scr,      &
                     exp_t_mat,  & 
                     ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  INTEGER                                  :: i
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(2,2)                   :: exp_t_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  DO i=1,ni
     v_scr(i,1) = v(i,1)  * exp_t_mat(1,1)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,1)    
 
     v_scr(i,2) = v(i,1)  * exp_t_mat(1,2)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,2)    

  END DO
!
!       Copy the temporary vector back to the input vector.
!
  v(:,1:2) = v_scr(:,1:2) 
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
  SUBROUTINE v_v_m_3(v,          &
                     v_scr,      &
                     exp_t_mat,  & 
                     ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  INTEGER                                  :: i
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(3,3)                   :: exp_t_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       
  DO i=1,ni
     v_scr(i,1) = v(i,1)  * exp_t_mat(1,1)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,1)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,1)    

     v_scr(i,2) = v(i,1)  * exp_t_mat(1,2)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,2)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,2)    

     v_scr(i,3) = v(i,1)  * exp_t_mat(1,3)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,3)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,3)    

  END DO    
!
!       Copy the temporary vector back to the input vector.
!
  v(:,1:3) = v_scr(:,1:3) 
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
  SUBROUTINE v_v_m_4(v,          &
                     v_scr,      &
                     exp_t_mat,  & 
                     ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  INTEGER                                  :: i
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(4,4)                   :: exp_t_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
  DO i=1,ni
     v_scr(i,1) = v(i,1)  * exp_t_mat(1,1)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,1)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,1)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,1)    

     v_scr(i,2) = v(i,1)  * exp_t_mat(1,2)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,2)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,2)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,2)    

     v_scr(i,3) = v(i,1)  * exp_t_mat(1,3)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,3)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,3)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,3)    

     v_scr(i,4) = v(i,1)  * exp_t_mat(1,4)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,4)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,4)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,4)    

  END DO
!
!       Copy the temporary vector back to the input vector.
!
  v(:,1:4) = v_scr(:,1:4) 
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
  SUBROUTINE v_v_m_5(v,          &
                     v_scr,      &
                     exp_t_mat,  & 
                     ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  INTEGER                                  :: i
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(5,5)                   :: exp_t_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  DO i=1,ni
     v_scr(i,1) = v(i,1)  * exp_t_mat(1,1)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,1)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,1)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,1)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,1)    

     v_scr(i,2) = v(i,1)  * exp_t_mat(1,2)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,2)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,2)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,2)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,2)    

     v_scr(i,3) = v(i,1)  * exp_t_mat(1,3)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,3)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,3)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,3)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,3)    

     v_scr(i,4) = v(i,1)  * exp_t_mat(1,4)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,4)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,4)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,4)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,4)    

     v_scr(i,5) = v(i,1)  * exp_t_mat(1,5)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,5)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,5)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,5)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,5)    

  END DO
!
!       Copy the temporary vector back to the input vector.
!
  v(:,1:5) = v_scr(:,1:5) 
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
  SUBROUTINE v_v_m_6(v,          &
                     v_scr,      &
                     exp_t_mat,  & 
                     ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  INTEGER                                  :: i
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(6,6)                   :: exp_t_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
  DO i=1,ni
     v_scr(i,1) = v(i,1)  * exp_t_mat(1,1)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,1)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,1)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,1)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,1)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,1)    

     v_scr(i,2) = v(i,1)  * exp_t_mat(1,2)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,2)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,2)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,2)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,2)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,2)    

     v_scr(i,3) = v(i,1)  * exp_t_mat(1,3)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,3)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,3)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,3)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,3)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,3)    

     v_scr(i,4) = v(i,1)  * exp_t_mat(1,4)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,4)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,4)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,4)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,4)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,4)    

     v_scr(i,5) = v(i,1)  * exp_t_mat(1,5)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,5)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,5)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,5)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,5)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,5)    

     v_scr(i,6) = v(i,1)  * exp_t_mat(1,6)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,6)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,6)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,6)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,6)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,6)    

  END DO                  
!
!       Copy the temporary vector back to the input vector.
!
  v(:,1:6) = v_scr(:,1:6) 
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
  SUBROUTINE v_v_m_7(v,          &
                     v_scr,      &
                     exp_t_mat,  & 
                     ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  INTEGER                                  :: i
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(7,7)                   :: exp_t_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  DO i=1,ni
     v_scr(i,1) = v(i,1)  * exp_t_mat(1,1)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,1)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,1)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,1)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,1)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,1)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,1)    

     v_scr(i,2) = v(i,1)  * exp_t_mat(1,2)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,2)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,2)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,2)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,2)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,2)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,2)    

     v_scr(i,3) = v(i,1)  * exp_t_mat(1,3)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,3)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,3)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,3)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,3)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,3)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,3)    

     v_scr(i,4) = v(i,1)  * exp_t_mat(1,4)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,4)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,4)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,4)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,4)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,4)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,4)    

     v_scr(i,5) = v(i,1)  * exp_t_mat(1,5)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,5)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,5)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,5)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,5)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,5)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,5)    

     v_scr(i,6) = v(i,1)  * exp_t_mat(1,6)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,6)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,6)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,6)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,6)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,6)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,6)    

     v_scr(i,7) = v(i,1)  * exp_t_mat(1,7)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,7)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,7)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,7)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,7)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,7)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,7)    

  END DO          
!
!       Copy the temporary vector back to the input vector.
!
  v(:,1:7) = v_scr(:,1:7) 
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
  SUBROUTINE v_v_m_8(v,          &
                     v_scr,      &
                     exp_t_mat,  & 
                     ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  INTEGER                                  :: i
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(8,8)                   :: exp_t_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  DO i=1,ni
     v_scr(i,1) = v(i,1)  * exp_t_mat(1,1)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,1)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,1)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,1)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,1)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,1)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,1)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,1)    

     v_scr(i,2) = v(i,1)  * exp_t_mat(1,2)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,2)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,2)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,2)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,2)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,2)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,2)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,2)    

     v_scr(i,3) = v(i,1)  * exp_t_mat(1,3)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,3)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,3)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,3)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,3)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,3)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,3)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,2)    

     v_scr(i,4) = v(i,1)  * exp_t_mat(1,4)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,4)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,4)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,4)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,4)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,4)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,4)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,4)    

     v_scr(i,5) = v(i,1)  * exp_t_mat(1,5)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,5)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,5)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,5)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,5)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,5)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,5)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,5)    

     v_scr(i,6) = v(i,1)  * exp_t_mat(1,6)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,6)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,6)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,6)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,6)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,6)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,6)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,6)    

     v_scr(i,7) = v(i,1)  * exp_t_mat(1,7)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,7)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,7)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,7)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,7)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,7)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,7)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,7)    

     v_scr(i,8) = v(i,1)  * exp_t_mat(1,8)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,8)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,8)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,8)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,8)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,8)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,8)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,8)    

  END DO
!
!       Copy the temporary vector back to the input vector.
!
  v(:,1:8) = v_scr(:,1:8) 
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
  SUBROUTINE v_v_m_9(v,          &
                     v_scr,      &
                     exp_t_mat,  & 
                     ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  INTEGER                                  :: i
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(9,9)                   :: exp_t_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  DO i=1,ni
     v_scr(i,1) = v(i,1)  * exp_t_mat(1,1)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,1)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,1)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,1)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,1)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,1)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,1)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,1)    &
                          +                   &
                  v(i,9)  * exp_t_mat(9,1)    

     v_scr(i,2) = v(i,1)  * exp_t_mat(1,2)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,2)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,2)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,2)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,2)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,2)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,2)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,2)    &
                          +                   &
                  v(i,9)  * exp_t_mat(9,2)    

     v_scr(i,3) = v(i,1)  * exp_t_mat(1,3)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,3)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,3)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,3)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,3)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,3)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,3)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,3)    &
                          +                   &
                  v(i,9)  * exp_t_mat(9,3)    

     v_scr(i,4) = v(i,1)  * exp_t_mat(1,4)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,4)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,4)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,4)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,4)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,4)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,4)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,4)    &
                          +                   &
                  v(i,9)  * exp_t_mat(9,4)    

     v_scr(i,5) = v(i,1)  * exp_t_mat(1,5)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,5)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,5)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,5)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,5)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,5)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,5)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,5)    &
                          +                   &
                  v(i,9)  * exp_t_mat(9,5)    

     v_scr(i,6) = v(i,1)  * exp_t_mat(1,6)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,6)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,6)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,6)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,6)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,6)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,6)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,6)    &
                          +                   &
                  v(i,9)  * exp_t_mat(9,6)    

     v_scr(i,7) = v(i,1)  * exp_t_mat(1,7)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,7)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,7)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,7)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,7)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,7)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,7)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,7)    &
                          +                   &
                  v(i,9)  * exp_t_mat(9,7)    

     v_scr(i,8) = v(i,1)  * exp_t_mat(1,8)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,8)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,8)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,8)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,8)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,8)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,8)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,8)    &
                          +                   &
                  v(i,9)  * exp_t_mat(9,8)    

     v_scr(i,9) = v(i,1)  * exp_t_mat(1,9)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,9)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,9)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,9)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,9)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,9)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,9)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,9)    &
                          +                   &
                  v(i,9)  * exp_t_mat(9,9)    
  END DO
!
!       Copy the temporary vector back to the input vector.
!
  v(:,1:9) = v_scr(:,1:9) 
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
  SUBROUTINE v_v_m_10(v,          &
                      v_scr,      &
                      exp_t_mat,  & 
                      ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  INTEGER                                  :: i
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(10,10)                 :: exp_t_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  DO i=1,ni
     v_scr(i,1) = v(i,1)  * exp_t_mat(1,1)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,1)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,1)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,1)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,1)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,1)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,1)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,1)    &
                          +                   &
                  v(i,9)  * exp_t_mat(9,1)    &
                          +                   &
                  v(i,10) * exp_t_mat(10,1)   

     v_scr(i,2) = v(i,1)  * exp_t_mat(1,2)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,2)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,2)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,2)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,2)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,2)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,2)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,2)    &
                          +                   &
                  v(i,9)  * exp_t_mat(9,2)    &
                          +                   &
                  v(i,10) * exp_t_mat(10,2)   

     v_scr(i,3) = v(i,1)  * exp_t_mat(1,3)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,3)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,3)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,3)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,3)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,3)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,3)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,3)    &
                          +                   &
                  v(i,9)  * exp_t_mat(9,3)    &
                          +                   &
                  v(i,10) * exp_t_mat(10,3)   

     v_scr(i,4) = v(i,1)  * exp_t_mat(1,4)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,4)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,4)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,4)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,4)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,4)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,4)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,4)    &
                          +                   &
                  v(i,9)  * exp_t_mat(9,4)    &
                          +                   &
                  v(i,10) * exp_t_mat(10,4)   

     v_scr(i,5) = v(i,1)  * exp_t_mat(1,5)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,5)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,5)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,5)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,5)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,5)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,5)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,5)    &
                          +                   &
                  v(i,9)  * exp_t_mat(9,5)    &
                          +                   &
                  v(i,10) * exp_t_mat(10,5)   

     v_scr(i,6) = v(i,1)  * exp_t_mat(1,6)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,6)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,6)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,6)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,6)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,6)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,6)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,6)    &
                          +                   &
                  v(i,9)  * exp_t_mat(9,6)    &
                          +                   &
                  v(i,10) * exp_t_mat(10,6)   

     v_scr(i,7) = v(i,1)  * exp_t_mat(1,7)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,7)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,7)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,7)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,7)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,7)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,7)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,7)    &
                          +                   &
                  v(i,9)  * exp_t_mat(9,7)    &
                          +                   &
                  v(i,10) * exp_t_mat(10,7)   

     v_scr(i,8) = v(i,1)  * exp_t_mat(1,8)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,8)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,8)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,8)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,8)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,8)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,8)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,8)    &
                          +                   &
                  v(i,9)  * exp_t_mat(9,8)    &
                          +                   &
                  v(i,10) * exp_t_mat(10,8)   

     v_scr(i,9) = v(i,1)  * exp_t_mat(1,9)    &
                          +                   &
                  v(i,2)  * exp_t_mat(2,9)    &
                          +                   &
                  v(i,3)  * exp_t_mat(3,9)    &
                          +                   &
                  v(i,4)  * exp_t_mat(4,9)    &
                          +                   &
                  v(i,5)  * exp_t_mat(5,9)    &
                          +                   &
                  v(i,6)  * exp_t_mat(6,9)    &
                          +                   &
                  v(i,7)  * exp_t_mat(7,9)    &
                          +                   &
                  v(i,8)  * exp_t_mat(8,9)    &
                          +                   &
                  v(i,9)  * exp_t_mat(9,9)    &
                          +                   &
                  v(i,10) * exp_t_mat(10,9)   

     v_scr(i,10) = v(i,1)  * exp_t_mat(1,10)    &
                          +                    &
                  v(i,2)  * exp_t_mat(2,10)    &
                          +                    &
                  v(i,3)  * exp_t_mat(3,10)    &
                          +                    &
                  v(i,4)  * exp_t_mat(4,10)    &
                          +                    &
                  v(i,5)  * exp_t_mat(5,10)    &
                          +                    &
                  v(i,6)  * exp_t_mat(6,10)    &
                          +                    &
                  v(i,7)  * exp_t_mat(7,10)    &
                          +                    &
                  v(i,8)  * exp_t_mat(8,10)    &
                          +                    &
                  v(i,9)  * exp_t_mat(9,10)    &
                          +                    &
                  v(i,10) * exp_t_mat(10,10)   
  END DO
!
!       Copy the temporary vector back to the input vector.
!
  v(:,1:10) = v_scr(:,1:10) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_m_10
!
END MODULE v_v_m
