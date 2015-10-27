!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         MODULE v_v_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      INTERFACE v_v_mat_gen
                 MODULE PROCEDURE v_v_mat_gen_d,                     &
                                  v_v_mat_gen_z
                      END INTERFACE  v_v_mat_gen
                      INTERFACE v_v_mat_2
                 MODULE PROCEDURE v_v_mat_2_d,                       &
                                  v_v_mat_2_z
                      END INTERFACE  v_v_mat_2
                      INTERFACE v_v_mat_3
                 MODULE PROCEDURE v_v_mat_3_d,                       &
                                  v_v_mat_3_z
                      END INTERFACE  v_v_mat_3
                      INTERFACE v_v_mat_4
                 MODULE PROCEDURE v_v_mat_4_d,                       &
                                  v_v_mat_4_z
                      END INTERFACE  v_v_mat_4
                      INTERFACE v_v_mat_5
                 MODULE PROCEDURE v_v_mat_5_d,                       &
                                  v_v_mat_5_z
                      END INTERFACE  v_v_mat_5
                      INTERFACE v_v_mat_6
                 MODULE PROCEDURE v_v_mat_6_d,                       &
                                  v_v_mat_6_z
                      END INTERFACE  v_v_mat_6
                      INTERFACE v_v_mat_7
                 MODULE PROCEDURE v_v_mat_7_d,                       &
                                  v_v_mat_7_z
                      END INTERFACE  v_v_mat_7
                      INTERFACE v_v_mat_8
                 MODULE PROCEDURE v_v_mat_8_d,                       &
                                  v_v_mat_8_z
                      END INTERFACE  v_v_mat_8
                      INTERFACE v_v_mat_9
                 MODULE PROCEDURE v_v_mat_9_d,                       &
                                  v_v_mat_9_z
                      END INTERFACE  v_v_mat_9
                      INTERFACE v_v_mat_10
                 MODULE PROCEDURE v_v_mat_10_d,                      &
                                  v_v_mat_10_z
                      END INTERFACE  v_v_mat_10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**begin prologue     v_v_mat
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       v_v_mat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                 CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck v_v_mat_gen_d
!***begin prologue     v_v_mat_gen_d     
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
  SUBROUTINE v_v_mat_gen_d(v,          &
                           v_scr,      &
                           dvr_mat,    & 
                           ni,nj,nk)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj, nk
  INTEGER                              :: j, k
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(nk,nk)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
!
!
  DO k=1,nk
     DO j=1,nk
        v_scr(:,j) = v_scr(:,j) + v(:,k) * dvr_mat(k,j) 
      END DO
   END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE v_v_mat_gen_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_v_mat_2_d
!***begin prologue     v_v_mat_2_d     
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
!***end prologue       v_v_mat_2_d
  SUBROUTINE v_v_mat_2_d(v,             &
                         v_scr,         &
                         dvr_mat,       & 
                         ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(2,2)               :: dvr_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  v_scr(:,1) = v_scr(:,1)                &
                       +                 &
               v(:,1)  * dvr_mat(1,1)    &
                       +                 &
               v(:,2)  * dvr_mat(2,1)    
 
  v_scr(:,2) = v_scr(:,2)                &
                       +                 &
               v(:,1)  * dvr_mat(1,2)    &
                       +                 &
               v(:,2)  * dvr_mat(2,2)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_2_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_v_mat_3_d
!***begin prologue     v_v_mat_3_d     
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
!***end prologue       v_v_mat_3_d
!
  SUBROUTINE v_v_mat_3_d(v,          &
                         v_scr,      &
                         dvr_mat,    & 
                         ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(3,3)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       
  v_scr(:,1) = v_scr(:,1)                &
                       +                 &
               v(:,1)  * dvr_mat(1,1)    &
                       +                 &
               v(:,2)  * dvr_mat(2,1)    &
                       +                 &
               v(:,3)  * dvr_mat(3,1)    

  v_scr(:,2) = v_scr(:,2)                &
                       +                 &
               v(:,1)  * dvr_mat(1,2)    &
                       +                 &
               v(:,2)  * dvr_mat(2,2)    &
                       +                 &
               v(:,3)  * dvr_mat(3,2)    

  v_scr(:,3) = v_scr(:,3)                &
                       +                 &
               v(:,1)  * dvr_mat(1,3)    &
                       +                 &
               v(:,2)  * dvr_mat(2,3)    &
                       +                 &
               v(:,3)  * dvr_mat(3,3)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_3_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_v_mat_4_d
!***begin prologue     v_v_mat_4_d     
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
!***end prologue       v_v_mat_4_d
!
  SUBROUTINE v_v_mat_4_d(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(4,4)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
  v_scr(:,1) = v_scr(:,1)                &
                       +                 &
               v(:,1)  * dvr_mat(1,1)    &
                       +                 &
               v(:,2)  * dvr_mat(2,1)    &
                       +                 &
               v(:,3)  * dvr_mat(3,1)    &
                       +                 &
               v(:,4)  * dvr_mat(4,1)    

  v_scr(:,2) = v_scr(:,2)                &
                       +                 &
               v(:,1)  * dvr_mat(1,2)    &
                       +                 &
               v(:,2)  * dvr_mat(2,2)    &
                       +                 &
               v(:,3)  * dvr_mat(3,2)    &
                       +                 &
               v(:,4)  * dvr_mat(4,2)    

  v_scr(:,3) = v_scr(:,3)                &
                       +                 &
               v(:,1)  * dvr_mat(1,3)    &
                       +                 &
               v(:,2)  * dvr_mat(2,3)    &
                       +                 &
               v(:,3)  * dvr_mat(3,3)    &
                       +                 &
               v(:,4)  * dvr_mat(4,3)    

  v_scr(:,4) = v_scr(:,4)                &
                       +                 &
               v(:,1)  * dvr_mat(1,4)    &
                       +                 &
               v(:,2)  * dvr_mat(2,4)    &
                       +                 &
               v(:,3)  * dvr_mat(3,4)    &
                       +                 &
               v(:,4)  * dvr_mat(4,4)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_4_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_v_mat_5_d
!***begin prologue     v_v_mat_5_d     
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
!***end prologue       v_v_mat_5_d
!
  SUBROUTINE v_v_mat_5_d(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(5,5)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_5_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_v_mat_6_d
!***begin prologue     v_v_mat_6_d     
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
!***end prologue       v_v_mat_6_d
!
  SUBROUTINE v_v_mat_6_d(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(6,6)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)    

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_6_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_v_mat_7_d
!***begin prologue     v_v_mat_7_d     
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
!***end prologue       v_v_mat_7_d
!
  SUBROUTINE v_v_mat_7_d(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(7,7)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)    

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)    

  v_scr(:,7) = v_scr(:,7)              &
                       +               &
               v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_7_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_v_mat_8_d
!***begin prologue     v_v_mat_8_d     
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
!***end prologue       v_v_mat_8_d
!
  SUBROUTINE v_v_mat_8_d(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(8,8)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)  &
                       +               &
               v(:,8)  * dvr_mat(8,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)  &
                       +               &
               v(:,8)  * dvr_mat(8,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)  &
                       +               &
               v(:,8)  * dvr_mat(8,5)    

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)  &
                       +               &
               v(:,8)  * dvr_mat(8,6)    

  v_scr(:,7) = v_scr(:,7)              &
                       +               &
               v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)  &
                       +               &
               v(:,8)  * dvr_mat(8,7)    

  v_scr(:,8) = v_scr(:,8)              &
                       +               &
               v(:,1)  * dvr_mat(1,8)  &
                       +               &
               v(:,2)  * dvr_mat(2,8)  &
                       +               &
               v(:,3)  * dvr_mat(3,8)  &
                       +               &
               v(:,4)  * dvr_mat(4,8)  &
                       +               &
               v(:,5)  * dvr_mat(5,8)  &
                       +               &
               v(:,6)  * dvr_mat(6,8)  &
                       +               &
               v(:,7)  * dvr_mat(7,8)  &
                       +               &
               v(:,8)  * dvr_mat(8,8)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_8_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_9_d
!***begin prologue     v_v_mat_9_d     
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
!***end prologue       v_v_mat_9_d
!
  SUBROUTINE v_v_mat_9_d(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(9,9)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)  &
                       +               &
               v(:,8)  * dvr_mat(8,1)  &
                       +               &
               v(:,9)  * dvr_mat(9,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)  &
                       +               &
               v(:,9)  * dvr_mat(9,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)  &
                       +               &
               v(:,8)  * dvr_mat(8,3)  &
                       +               &
               v(:,9)  * dvr_mat(9,3)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)  &
                       +               &
               v(:,8)  * dvr_mat(8,4)  &
                       +               &
               v(:,9)  * dvr_mat(9,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)  &
                       +               &
               v(:,8)  * dvr_mat(8,5)  &
                       +               &
               v(:,9)  * dvr_mat(9,5)    

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)  &
                       +               &
               v(:,8)  * dvr_mat(8,6)  &
                       +               &
               v(:,9)  * dvr_mat(9,6)    

  v_scr(:,7) = v_scr(:,7)              &
                       +               &
               v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)  &
                       +               &
               v(:,8)  * dvr_mat(8,7)  &
                       +               &
               v(:,9)  * dvr_mat(9,7)    

  v_scr(:,8) = v_scr(:,8)              &
                       +               &
               v(:,1)  * dvr_mat(1,8)  &
                       +               &
               v(:,2)  * dvr_mat(2,8)  &
                       +               &
               v(:,3)  * dvr_mat(3,8)  &
                       +               &
               v(:,4)  * dvr_mat(4,8)  &
                       +               &
               v(:,5)  * dvr_mat(5,8)  &
                       +               &
               v(:,6)  * dvr_mat(6,8)  &
                       +               &
               v(:,7)  * dvr_mat(7,8)  &
                       +               &
               v(:,8)  * dvr_mat(8,8)  &
                       +               &
               v(:,9)  * dvr_mat(9,8)    

  v_scr(:,9) = v_scr(:,9)              &
                       +               &
               v(:,1)  * dvr_mat(1,9)  &
                       +               &
               v(:,2)  * dvr_mat(2,9)  &
                       +               &
               v(:,3)  * dvr_mat(3,9)  &
                       +               &
               v(:,4)  * dvr_mat(4,9)  &
                       +               &
               v(:,5)  * dvr_mat(5,9)  &
                       +               &
               v(:,6)  * dvr_mat(6,9)  &
                       +               &
               v(:,7)  * dvr_mat(7,9)  &
                       +               &
               v(:,8)  * dvr_mat(8,9)  &
                       +               &
               v(:,9)  * dvr_mat(9,9)    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_9_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_v_mat_10_d
!***begin prologue     v_v_mat_10_d     
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
!***end prologue       v_v_mat_10_d
!
  SUBROUTINE v_v_mat_10_d(v,        &
                          v_scr,    &
                          dvr_mat,  & 
                          ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(10,10)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)  &
                       +               &
               v(:,8)  * dvr_mat(8,1)  &
                       +               &
               v(:,9)  * dvr_mat(9,1)  &
                       +               &
               v(:,10) * dvr_mat(10,1)   

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)  &
                       +               &
               v(:,9)  * dvr_mat(9,2)  &
                       +               &
               v(:,10) * dvr_mat(10,2)   

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)  &
                       +               &
               v(:,8)  * dvr_mat(8,3)  &
                       +               &
               v(:,9)  * dvr_mat(9,3)  &
                       +               &
               v(:,10) * dvr_mat(10,3)   

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)  &
                       +               &
               v(:,8)  * dvr_mat(8,4)  &
                       +               &
               v(:,9)  * dvr_mat(9,4)  &
                       +               &
               v(:,10) * dvr_mat(10,4)   

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)  &
                       +               &
               v(:,8)  * dvr_mat(8,5)  &
                       +               &
               v(:,9)  * dvr_mat(9,5)  &
                       +               &
               v(:,10) * dvr_mat(10,5)   

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)  &
                       +               &
               v(:,8)  * dvr_mat(8,6)  &
                       +               &
               v(:,9)  * dvr_mat(9,6)  &
                       +               &
               v(:,10) * dvr_mat(10,6)   

  v_scr(:,7) = v_scr(:,7)              &
                       +               &
               v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)  &
                       +               &
               v(:,8)  * dvr_mat(8,7)  &
                       +               &
               v(:,9)  * dvr_mat(9,7)  &
                       +               &
               v(:,10) * dvr_mat(10,7)   

  v_scr(:,8) = v_scr(:,8)              &
                       +               &
               v(:,1)  * dvr_mat(1,8)  &
                       +               &
               v(:,2)  * dvr_mat(2,8)  &
                       +               &
               v(:,3)  * dvr_mat(3,8)  &
                       +               &
               v(:,4)  * dvr_mat(4,8)  &
                       +               &
               v(:,5)  * dvr_mat(5,8)  &
                       +               &
               v(:,6)  * dvr_mat(6,8)  &
                       +               &
               v(:,7)  * dvr_mat(7,8)  &
                       +               &
               v(:,8)  * dvr_mat(8,8)  &
                       +               &
               v(:,9)  * dvr_mat(9,8)  &
                       +               &
               v(:,10) * dvr_mat(10,8)   

  v_scr(:,9) = v_scr(:,9)              &
                       +               &
               v(:,1)  * dvr_mat(1,9)  &
                       +               &
               v(:,2)  * dvr_mat(2,9)  &
                       +               &
               v(:,3)  * dvr_mat(3,9)  &
                       +               &
               v(:,4)  * dvr_mat(4,9)  &
                       +               &
               v(:,5)  * dvr_mat(5,9)  &
                       +               &
               v(:,6)  * dvr_mat(6,9)  &
                       +               &
               v(:,7)  * dvr_mat(7,9)  &
                       +               &
               v(:,8)  * dvr_mat(8,9)  &
                       +               &
               v(:,9)  * dvr_mat(9,9)  &
                       +               &
               v(:,10) * dvr_mat(10,9)   

  v_scr(:,10) = v_scr(:,10)            &
                       +               &
               v(:,1)  * dvr_mat(1,10)  &
                       +                &
               v(:,2)  * dvr_mat(2,10)  &
                       +                &
               v(:,3)  * dvr_mat(3,10)  &
                       +                &
               v(:,4)  * dvr_mat(4,10)  &
                       +                &
               v(:,5)  * dvr_mat(5,10)  &
                       +                &
               v(:,6)  * dvr_mat(6,10)  &
                       +                &
               v(:,7)  * dvr_mat(7,10)  &
                       +                &
               v(:,8)  * dvr_mat(8,10)  &
                       +                &
               v(:,9)  * dvr_mat(9,10)  &
                       +                &
               v(:,10) * dvr_mat(10,10)   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_10_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck v_v_mat_gen_z
!***begin prologue     v_v_mat_gen_z     
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
  SUBROUTINE v_v_mat_gen_z(v,          &
                           v_scr,      &
                           dvr_mat,    & 
                           ni,nj,nk)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj, nk
  INTEGER                              :: j, k
  COMPLEX*16, DIMENSION(:,:)           :: v
  COMPLEX*16, DIMENSION(:,:)           :: v_scr
  REAL*8, DIMENSION(nk,nk)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
!
!
  DO k=1,nk
     DO j=1,nk
        v_scr(:,j) = v_scr(:,j) + v(:,k) * dvr_mat(k,j) 
      END DO
   END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE v_v_mat_gen_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_v_mat_2_z
!***begin prologue     v_v_mat_2_z     
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
!***end prologue       v_v_mat_2_z
  SUBROUTINE v_v_mat_2_z(v,             &
                         v_scr,         &
                         dvr_mat,       & 
                         ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(:,:)           :: v
  COMPLEX*16, DIMENSION(:,:)           :: v_scr
  REAL*8, DIMENSION(2,2)               :: dvr_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  v_scr(:,1) = v_scr(:,1)                &
                       +                 &
               v(:,1)  * dvr_mat(1,1)    &
                       +                 &
               v(:,2)  * dvr_mat(2,1)    
 
  v_scr(:,2) = v_scr(:,2)                &
                       +                 &
               v(:,1)  * dvr_mat(1,2)    &
                       +                 &
               v(:,2)  * dvr_mat(2,2)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_2_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_v_mat_3_z
!***begin prologue     v_v_mat_3_z     
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
!***end prologue       v_v_mat_3_z
!
  SUBROUTINE v_v_mat_3_z(v,          &
                         v_scr,      &
                         dvr_mat,    & 
                         ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(:,:)           :: v
  COMPLEX*16, DIMENSION(:,:)           :: v_scr
  REAL*8, DIMENSION(3,3)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       
  v_scr(:,1) = v_scr(:,1)                &
                       +                 &
               v(:,1)  * dvr_mat(1,1)    &
                       +                 &
               v(:,2)  * dvr_mat(2,1)    &
                       +                 &
               v(:,3)  * dvr_mat(3,1)    

  v_scr(:,2) = v_scr(:,2)                &
                       +                 &
               v(:,1)  * dvr_mat(1,2)    &
                       +                 &
               v(:,2)  * dvr_mat(2,2)    &
                       +                 &
               v(:,3)  * dvr_mat(3,2)    

  v_scr(:,3) = v_scr(:,3)                &
                       +                 &
               v(:,1)  * dvr_mat(1,3)    &
                       +                 &
               v(:,2)  * dvr_mat(2,3)    &
                       +                 &
               v(:,3)  * dvr_mat(3,3)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_3_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_v_mat_4_z
!***begin prologue     v_v_mat_4_z     
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
!***end prologue       v_v_mat_4_z
!
  SUBROUTINE v_v_mat_4_z(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(:,:)           :: v
  COMPLEX*16, DIMENSION(:,:)           :: v_scr
  REAL*8, DIMENSION(4,4)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
  v_scr(:,1) = v_scr(:,1)                &
                       +                 &
               v(:,1)  * dvr_mat(1,1)    &
                       +                 &
               v(:,2)  * dvr_mat(2,1)    &
                       +                 &
               v(:,3)  * dvr_mat(3,1)    &
                       +                 &
               v(:,4)  * dvr_mat(4,1)    

  v_scr(:,2) = v_scr(:,2)                &
                       +                 &
               v(:,1)  * dvr_mat(1,2)    &
                       +                 &
               v(:,2)  * dvr_mat(2,2)    &
                       +                 &
               v(:,3)  * dvr_mat(3,2)    &
                       +                 &
               v(:,4)  * dvr_mat(4,2)    

  v_scr(:,3) = v_scr(:,3)                &
                       +                 &
               v(:,1)  * dvr_mat(1,3)    &
                       +                 &
               v(:,2)  * dvr_mat(2,3)    &
                       +                 &
               v(:,3)  * dvr_mat(3,3)    &
                       +                 &
               v(:,4)  * dvr_mat(4,3)    

  v_scr(:,4) = v_scr(:,4)                &
                       +                 &
               v(:,1)  * dvr_mat(1,4)    &
                       +                 &
               v(:,2)  * dvr_mat(2,4)    &
                       +                 &
               v(:,3)  * dvr_mat(3,4)    &
                       +                 &
               v(:,4)  * dvr_mat(4,4)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_4_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_v_mat_5_z
!***begin prologue     v_v_mat_5_z     
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
!***end prologue       v_v_mat_5_z
!
  SUBROUTINE v_v_mat_5_z(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(:,:)           :: v
  COMPLEX*16, DIMENSION(:,:)           :: v_scr
  REAL*8, DIMENSION(5,5)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_5_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_v_mat_6_z
!***begin prologue     v_v_mat_6_z     
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
!***end prologue       v_v_mat_6_z
!
  SUBROUTINE v_v_mat_6_z(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(:,:)           :: v
  COMPLEX*16, DIMENSION(:,:)           :: v_scr
  REAL*8, DIMENSION(6,6)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)    

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_6_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_7_z
!***begin prologue     v_v_mat_7_z     
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
!***end prologue       v_v_mat_7_z
!
  SUBROUTINE v_v_mat_7_z(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(:,:)           :: v
  COMPLEX*16, DIMENSION(:,:)           :: v_scr
  REAL*8, DIMENSION(7,7)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)    

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)    

  v_scr(:,7) = v_scr(:,7)              &
                       +               &
               v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_7_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_v_mat_8_z
!***begin prologue     v_v_mat_8_z     
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
!***end prologue       v_v_mat_8_z
!
  SUBROUTINE v_v_mat_8_z(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(:,:)           :: v
  COMPLEX*16, DIMENSION(:,:)           :: v_scr
  REAL*8, DIMENSION(8,8)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)  &
                       +               &
               v(:,8)  * dvr_mat(8,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)  &
                       +               &
               v(:,8)  * dvr_mat(8,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)  &
                       +               &
               v(:,8)  * dvr_mat(8,5)    

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)  &
                       +               &
               v(:,8)  * dvr_mat(8,6)    

  v_scr(:,7) = v_scr(:,7)              &
                       +               &
               v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)  &
                       +               &
               v(:,8)  * dvr_mat(8,7)    

  v_scr(:,8) = v_scr(:,8)              &
                       +               &
               v(:,1)  * dvr_mat(1,8)  &
                       +               &
               v(:,2)  * dvr_mat(2,8)  &
                       +               &
               v(:,3)  * dvr_mat(3,8)  &
                       +               &
               v(:,4)  * dvr_mat(4,8)  &
                       +               &
               v(:,5)  * dvr_mat(5,8)  &
                       +               &
               v(:,6)  * dvr_mat(6,8)  &
                       +               &
               v(:,7)  * dvr_mat(7,8)  &
                       +               &
               v(:,8)  * dvr_mat(8,8)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_8_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_v_mat_9_z
!***begin prologue     v_v_mat_9_z     
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
!***end prologue       v_v_mat_9_z
!
  SUBROUTINE v_v_mat_9_z(v,        &
                         v_scr,    &
                         dvr_mat,  & 
                         ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(:,:)           :: v
  COMPLEX*16, DIMENSION(:,:)           :: v_scr
  REAL*8, DIMENSION(9,9)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)  &
                       +               &
               v(:,8)  * dvr_mat(8,1)  &
                       +               &
               v(:,9)  * dvr_mat(9,1)    

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)  &
                       +               &
               v(:,9)  * dvr_mat(9,2)    

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)  &
                       +               &
               v(:,8)  * dvr_mat(8,3)  &
                       +               &
               v(:,9)  * dvr_mat(9,3)    

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)  &
                       +               &
               v(:,8)  * dvr_mat(8,4)  &
                       +               &
               v(:,9)  * dvr_mat(9,4)    

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)  &
                       +               &
               v(:,8)  * dvr_mat(8,5)  &
                       +               &
               v(:,9)  * dvr_mat(9,5)    

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)  &
                       +               &
               v(:,8)  * dvr_mat(8,6)  &
                       +               &
               v(:,9)  * dvr_mat(9,6)    

  v_scr(:,7) = v_scr(:,7)              &
                       +               &
               v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)  &
                       +               &
               v(:,8)  * dvr_mat(8,7)  &
                       +               &
               v(:,9)  * dvr_mat(9,7)    

  v_scr(:,8) = v_scr(:,8)              &
                       +               &
               v(:,1)  * dvr_mat(1,8)  &
                       +               &
               v(:,2)  * dvr_mat(2,8)  &
                       +               &
               v(:,3)  * dvr_mat(3,8)  &
                       +               &
               v(:,4)  * dvr_mat(4,8)  &
                       +               &
               v(:,5)  * dvr_mat(5,8)  &
                       +               &
               v(:,6)  * dvr_mat(6,8)  &
                       +               &
               v(:,7)  * dvr_mat(7,8)  &
                       +               &
               v(:,8)  * dvr_mat(8,8)  &
                       +               &
               v(:,9)  * dvr_mat(9,8)    

  v_scr(:,9) = v_scr(:,9)              &
                       +               &
               v(:,1)  * dvr_mat(1,9)  &
                       +               &
               v(:,2)  * dvr_mat(2,9)  &
                       +               &
               v(:,3)  * dvr_mat(3,9)  &
                       +               &
               v(:,4)  * dvr_mat(4,9)  &
                       +               &
               v(:,5)  * dvr_mat(5,9)  &
                       +               &
               v(:,6)  * dvr_mat(6,9)  &
                       +               &
               v(:,7)  * dvr_mat(7,9)  &
                       +               &
               v(:,8)  * dvr_mat(8,9)  &
                       +               &
               v(:,9)  * dvr_mat(9,9)    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_9_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_v_mat_10_z
!***begin prologue     v_v_mat_10_z     
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
!***end prologue       v_v_mat_10_z
!
  SUBROUTINE v_v_mat_10_z(v,        &
                          v_scr,    &
                          dvr_mat,  & 
                          ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(:,:)           :: v
  COMPLEX*16, DIMENSION(:,:)           :: v_scr
  REAL*8, DIMENSION(10,10)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(:,1) = v_scr(:,1)              &
                       +               &
               v(:,1)  * dvr_mat(1,1)  &
                       +               &
               v(:,2)  * dvr_mat(2,1)  &
                       +               &
               v(:,3)  * dvr_mat(3,1)  &
                       +               &
               v(:,4)  * dvr_mat(4,1)  &
                       +               &
               v(:,5)  * dvr_mat(5,1)  &
                       +               &
               v(:,6)  * dvr_mat(6,1)  &
                       +               &
               v(:,7)  * dvr_mat(7,1)  &
                       +               &
               v(:,8)  * dvr_mat(8,1)  &
                       +               &
               v(:,9)  * dvr_mat(9,1)  &
                       +               &
               v(:,10) * dvr_mat(10,1)   

  v_scr(:,2) = v_scr(:,2)              &
                       +               &
               v(:,1)  * dvr_mat(1,2)  &
                       +               &
               v(:,2)  * dvr_mat(2,2)  &
                       +               &
               v(:,3)  * dvr_mat(3,2)  &
                       +               &
               v(:,4)  * dvr_mat(4,2)  &
                       +               &
               v(:,5)  * dvr_mat(5,2)  &
                       +               &
               v(:,6)  * dvr_mat(6,2)  &
                       +               &
               v(:,7)  * dvr_mat(7,2)  &
                       +               &
               v(:,8)  * dvr_mat(8,2)  &
                       +               &
               v(:,9)  * dvr_mat(9,2)  &
                       +               &
               v(:,10) * dvr_mat(10,2)   

  v_scr(:,3) = v_scr(:,3)              &
                       +               &
               v(:,1)  * dvr_mat(1,3)  &
                       +               &
               v(:,2)  * dvr_mat(2,3)  &
                       +               &
               v(:,3)  * dvr_mat(3,3)  &
                       +               &
               v(:,4)  * dvr_mat(4,3)  &
                       +               &
               v(:,5)  * dvr_mat(5,3)  &
                       +               &
               v(:,6)  * dvr_mat(6,3)  &
                       +               &
               v(:,7)  * dvr_mat(7,3)  &
                       +               &
               v(:,8)  * dvr_mat(8,3)  &
                       +               &
               v(:,9)  * dvr_mat(9,3)  &
                       +               &
               v(:,10) * dvr_mat(10,3)   

  v_scr(:,4) = v_scr(:,4)              &
                       +               &
               v(:,1)  * dvr_mat(1,4)  &
                       +               &
               v(:,2)  * dvr_mat(2,4)  &
                       +               &
               v(:,3)  * dvr_mat(3,4)  &
                       +               &
               v(:,4)  * dvr_mat(4,4)  &
                       +               &
               v(:,5)  * dvr_mat(5,4)  &
                       +               &
               v(:,6)  * dvr_mat(6,4)  &
                       +               &
               v(:,7)  * dvr_mat(7,4)  &
                       +               &
               v(:,8)  * dvr_mat(8,4)  &
                       +               &
               v(:,9)  * dvr_mat(9,4)  &
                       +               &
               v(:,10) * dvr_mat(10,4)   

  v_scr(:,5) = v_scr(:,5)              &
                       +               &
               v(:,1)  * dvr_mat(1,5)  &
                       +               &
               v(:,2)  * dvr_mat(2,5)  &
                       +               &
               v(:,3)  * dvr_mat(3,5)  &
                       +               &
               v(:,4)  * dvr_mat(4,5)  &
                       +               &
               v(:,5)  * dvr_mat(5,5)  &
                       +               &
               v(:,6)  * dvr_mat(6,5)  &
                       +               &
               v(:,7)  * dvr_mat(7,5)  &
                       +               &
               v(:,8)  * dvr_mat(8,5)  &
                       +               &
               v(:,9)  * dvr_mat(9,5)  &
                       +               &
               v(:,10) * dvr_mat(10,5)   

  v_scr(:,6) = v_scr(:,6)              &
                       +               &
               v(:,1)  * dvr_mat(1,6)  &
                       +               &
               v(:,2)  * dvr_mat(2,6)  &
                       +               &
               v(:,3)  * dvr_mat(3,6)  &
                       +               &
               v(:,4)  * dvr_mat(4,6)  &
                       +               &
               v(:,5)  * dvr_mat(5,6)  &
                       +               &
               v(:,6)  * dvr_mat(6,6)  &
                       +               &
               v(:,7)  * dvr_mat(7,6)  &
                       +               &
               v(:,8)  * dvr_mat(8,6)  &
                       +               &
               v(:,9)  * dvr_mat(9,6)  &
                       +               &
               v(:,10) * dvr_mat(10,6)   

  v_scr(:,7) = v_scr(:,7)              &
                       +               &
               v(:,1)  * dvr_mat(1,7)  &
                       +               &
               v(:,2)  * dvr_mat(2,7)  &
                       +               &
               v(:,3)  * dvr_mat(3,7)  &
                       +               &
               v(:,4)  * dvr_mat(4,7)  &
                       +               &
               v(:,5)  * dvr_mat(5,7)  &
                       +               &
               v(:,6)  * dvr_mat(6,7)  &
                       +               &
               v(:,7)  * dvr_mat(7,7)  &
                       +               &
               v(:,8)  * dvr_mat(8,7)  &
                       +               &
               v(:,9)  * dvr_mat(9,7)  &
                       +               &
               v(:,10) * dvr_mat(10,7)   

  v_scr(:,8) = v_scr(:,8)              &
                       +               &
               v(:,1)  * dvr_mat(1,8)  &
                       +               &
               v(:,2)  * dvr_mat(2,8)  &
                       +               &
               v(:,3)  * dvr_mat(3,8)  &
                       +               &
               v(:,4)  * dvr_mat(4,8)  &
                       +               &
               v(:,5)  * dvr_mat(5,8)  &
                       +               &
               v(:,6)  * dvr_mat(6,8)  &
                       +               &
               v(:,7)  * dvr_mat(7,8)  &
                       +               &
               v(:,8)  * dvr_mat(8,8)  &
                       +               &
               v(:,9)  * dvr_mat(9,8)  &
                       +               &
               v(:,10) * dvr_mat(10,8)   

  v_scr(:,9) = v_scr(:,9)              &
                       +               &
               v(:,1)  * dvr_mat(1,9)  &
                       +               &
               v(:,2)  * dvr_mat(2,9)  &
                       +               &
               v(:,3)  * dvr_mat(3,9)  &
                       +               &
               v(:,4)  * dvr_mat(4,9)  &
                       +               &
               v(:,5)  * dvr_mat(5,9)  &
                       +               &
               v(:,6)  * dvr_mat(6,9)  &
                       +               &
               v(:,7)  * dvr_mat(7,9)  &
                       +               &
               v(:,8)  * dvr_mat(8,9)  &
                       +               &
               v(:,9)  * dvr_mat(9,9)  &
                       +               &
               v(:,10) * dvr_mat(10,9)   

  v_scr(:,10) = v_scr(:,10)            &
                       +               &
               v(:,1)  * dvr_mat(1,10)  &
                       +                &
               v(:,2)  * dvr_mat(2,10)  &
                       +                &
               v(:,3)  * dvr_mat(3,10)  &
                       +                &
               v(:,4)  * dvr_mat(4,10)  &
                       +                &
               v(:,5)  * dvr_mat(5,10)  &
                       +                &
               v(:,6)  * dvr_mat(6,10)  &
                       +                &
               v(:,7)  * dvr_mat(7,10)  &
                       +                &
               v(:,8)  * dvr_mat(8,10)  &
                       +                &
               v(:,9)  * dvr_mat(9,10)  &
                       +                &
               v(:,10) * dvr_mat(10,10)   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_v_mat_10_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                END      MODULE v_v_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
