!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         MODULE v_out_so_mat_v_in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      INTERFACE v_so_mat_v_gen
                 MODULE PROCEDURE v_so_mat_v_gen_d,                     &
                                  v_so_mat_v_gen_z
                      END INTERFACE  v_so_mat_v_gen
                      INTERFACE v_so_mat_v_2
                 MODULE PROCEDURE v_so_mat_v_2_d,                       &
                                  v_so_mat_v_2_z
                      END INTERFACE  v_so_mat_v_2
                      INTERFACE v_so_mat_v_3
                 MODULE PROCEDURE v_so_mat_v_3_d,                       &
                                  v_so_mat_v_3_z
                      END INTERFACE  v_so_mat_v_3
                      INTERFACE v_so_mat_v_4
                 MODULE PROCEDURE v_so_mat_v_4_d,                       &
                                  v_so_mat_v_4_z
                      END INTERFACE  v_so_mat_v_4
                      INTERFACE v_so_mat_v_5
                 MODULE PROCEDURE v_so_mat_v_5_d,                       &
                                  v_so_mat_v_5_z
                      END INTERFACE  v_so_mat_v_5
                      INTERFACE v_so_mat_v_6
                 MODULE PROCEDURE v_so_mat_v_6_d,                       &
                                  v_so_mat_v_6_z
                      END INTERFACE  v_so_mat_v_6
                      INTERFACE v_so_mat_v_7
                 MODULE PROCEDURE v_so_mat_v_7_d,                       &
                                  v_so_mat_v_7_z
                      END INTERFACE  v_so_mat_v_7
                      INTERFACE v_so_mat_v_8
                 MODULE PROCEDURE v_so_mat_v_8_d,                       &
                                  v_so_mat_v_8_z
                      END INTERFACE  v_so_mat_v_8
                      INTERFACE v_so_mat_v_9
                 MODULE PROCEDURE v_so_mat_v_9_d,                       &
                                  v_so_mat_v_9_z
                      END INTERFACE  v_so_mat_v_9
                      INTERFACE v_so_mat_v_10
                 MODULE PROCEDURE v_so_mat_v_10_d,                      &
                                  v_so_mat_v_10_z
                      END INTERFACE  v_so_mat_v_10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**begin prologue     v_out_so_mat_v_in
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       v_out_so_mat_v_in
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_gen_d
!***begin prologue     v_so_mat_v_gen_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propator matrix vector multiplies
!***                   for a general FEDVR Hamiltonian.
!***                   The parameter nj is a dummy and can take on values
!***                   consistent with any dimensional problem.  In a 1D
!***                   case nj=1, in 2D nj=nx and in 3D nj=ny*nx.
!
!***references
!***routines called    ebcxx, ambcxx, apbcxx
!***end prologue       
!
  SUBROUTINE v_so_mat_v_gen_d(v,             &
                              v_scr,         &
                              dvr_mat,       &
                              ni,nj,nk)
  IMPLICIT NONE
  INTEGER                              :: ni, nj, nk
  INTEGER                              :: i, k
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(nk,nk)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   DO i=1,nk
      DO k=1,nk
         v_scr(k,:) = v_scr(k,:) + dvr_mat(k,i) * v(i,:) 
      END DO
   END DO
!
   v(1:nk,:) = v_scr(1:nk,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE v_so_mat_v_gen_d
!
!deck v_so_mat_v_2_d
!***begin prologue     v_so_mat_v_2_d     
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
!***end prologue       v_so_mat_v_2_d
!
  SUBROUTINE v_so_mat_v_2_d(v,         &
                            v_scr,     &
                            dvr_mat,   &
                            ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(2,2)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   
  v_scr(2,:) = dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   
!
  v(1:2,:) = v_scr(1:2,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_2_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_so_mat_v_3_d
!***begin prologue     v_so_mat_v_3_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_3_d
!
  SUBROUTINE v_so_mat_v_3_d(v,            &
                            v_scr,        &
                            dvr_mat,      &
                            ni,nj)
  IMPLICIT NONE
  INTEGER                                :: ni, nj
  REAL*8, DIMENSION(:,:)                 :: v
  REAL*8, DIMENSION(:,:)                 :: v_scr
  REAL*8, DIMENSION(3,3)                 :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   
  v_scr(2,:) = dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   
  v_scr(3,:) = dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   
!
  v(1:3,:) = v_scr(1:3,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_3_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_so_mat_v_4_d
!***begin prologue     v_so_mat_v_4_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 4*4 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_4_d
!
  SUBROUTINE v_so_mat_v_4_d(v,             &
                            v_scr,         &
                            dvr_mat,       &
                            ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(4,4)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   &
                            +          &   
               dvr_mat(1,4) * v(4,:)   

  v_scr(2,:) = dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   &
                            +          &
               dvr_mat(2,4) * v(4,:)   

  v_scr(3,:) = dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   &
                            +          &
               dvr_mat(3,4) * v(4,:)   

  v_scr(4,:) = dvr_mat(4,1) * v(1,:)   &
                            +          &
               dvr_mat(4,2) * v(2,:)   &
                            +          &
               dvr_mat(4,3) * v(3,:)   &
                            +          &
               dvr_mat(4,4) * v(4,:)   
!
  v(1:4,:) = v_scr(1:4,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_4_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_so_mat_v_5_d
!***begin prologue     v_so_mat_v_5_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_5_d
!
  SUBROUTINE v_so_mat_v_5_d(v,          &
                            v_scr,      &
                            dvr_mat,    &
                            ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(5,5)               :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   &
                            +          &   
               dvr_mat(1,4) * v(4,:)   &
                            +          &
               dvr_mat(1,5) * v(5,:)      
  v_scr(2,:) = dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   &
                            +          &
               dvr_mat(2,4) * v(4,:)   &
                            +          &
               dvr_mat(2,5) * v(5,:)   
  v_scr(3,:) = dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   &
                            +          &
               dvr_mat(3,4) * v(4,:)   &
                            +          &
               dvr_mat(3,5) * v(5,:)   
  v_scr(4,:) = dvr_mat(4,1) * v(1,:)   &
                            +          &
               dvr_mat(4,2) * v(2,:)   &
                            +          &
               dvr_mat(4,3) * v(3,:)   &
                            +          &
               dvr_mat(4,4) * v(4,:)   &
                            +          &
               dvr_mat(4,5) * v(5,:)    
  v_scr(5,:) = dvr_mat(5,1) * v(1,:)   &
                            +          &
               dvr_mat(5,2) * v(2,:)   &
                            +          &
               dvr_mat(5,3) * v(3,:)   &
                            +          &
               dvr_mat(5,4) * v(4,:)   &
                            +          &
               dvr_mat(5,5) * v(5,:)    
!
  v(1:5,:) = v_scr(1:5,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_5_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_so_mat_v_6_d
!***begin prologue     v_so_mat_v_6_d     
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
!***end prologue       v_so_mat_v_6_d
!
  SUBROUTINE v_so_mat_v_6_d(v,           &    
                            v_scr,       &
                            dvr_mat,     &
                            ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(6,6)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   &
                            +          &   
               dvr_mat(1,4) * v(4,:)   &
                            +          &
               dvr_mat(1,5) * v(5,:)   &   
                            +          &
               dvr_mat(1,6) * v(6,:)      
  v_scr(2,:) = dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   &
                            +          &   
               dvr_mat(2,4) * v(4,:)   &
                            +          &
               dvr_mat(2,5) * v(5,:)   &   
                            +          &
               dvr_mat(2,6) * v(6,:)      
  v_scr(3,:) = dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   &
                            +          &   
               dvr_mat(3,4) * v(4,:)   &
                            +          &
               dvr_mat(3,5) * v(5,:)   &   
                            +          &
               dvr_mat(3,6) * v(6,:)      
  v_scr(4,:) = dvr_mat(4,1) * v(1,:)   &
                            +          &
               dvr_mat(4,2) * v(2,:)   &
                            +          &
               dvr_mat(4,3) * v(3,:)   &
                            +          &   
               dvr_mat(4,4) * v(4,:)   &
                            +          &
               dvr_mat(4,5) * v(5,:)   &   
                            +          &
               dvr_mat(4,6) * v(6,:)      
  v_scr(5,:) = dvr_mat(5,1) * v(1,:)   &
                            +          &
               dvr_mat(5,2) * v(2,:)   &
                            +          &
               dvr_mat(5,3) * v(3,:)   &
                            +          &   
               dvr_mat(5,4) * v(4,:)   &
                            +          &
               dvr_mat(5,5) * v(5,:)   &   
                            +          &
               dvr_mat(5,6) * v(6,:)      

  v_scr(6,:) = dvr_mat(6,1) * v(1,:)   &
                            +          &
               dvr_mat(6,2) * v(2,:)   &
                            +          &
               dvr_mat(6,3) * v(3,:)   &
                            +          &   
               dvr_mat(6,4) * v(4,:)   &
                            +          &
               dvr_mat(6,5) * v(5,:)   &   
                            +          &
               dvr_mat(6,6) * v(6,:)      
!
  v(1:6,:) = v_scr(1:6,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_6_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_so_mat_v_7_d
!***begin prologue     v_so_mat_v_7_d     
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
!***end prologue       v_so_mat_v_7_d
!
  SUBROUTINE v_so_mat_v_7_d(v,           &    
                            v_scr,       &
                            dvr_mat,     &
                            ni,nj)
                     
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(7,7)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)   &
                            +          &
               dvr_mat(1,2) * v(2,:)   &
                            +          &
               dvr_mat(1,3) * v(3,:)   &
                            +          &   
               dvr_mat(1,4) * v(4,:)   &
                            +          &
               dvr_mat(1,5) * v(5,:)   &   
                            +          &
               dvr_mat(1,6) * v(6,:)   &   
                            +          &
               dvr_mat(1,7) * v(7,:)      

  v_scr(2,:) = dvr_mat(2,1) * v(1,:)   &
                            +          &
               dvr_mat(2,2) * v(2,:)   &
                            +          &
               dvr_mat(2,3) * v(3,:)   &
                            +          &   
               dvr_mat(2,4) * v(4,:)   &
                            +          &
               dvr_mat(2,5) * v(5,:)   &   
                            +          &
               dvr_mat(2,6) * v(6,:)   &   
                            +          &
               dvr_mat(2,7) * v(7,:)      
  v_scr(3,:) = dvr_mat(3,1) * v(1,:)   &
                            +          &
               dvr_mat(3,2) * v(2,:)   &
                            +          &
               dvr_mat(3,3) * v(3,:)   &
                            +          &   
               dvr_mat(3,4) * v(4,:)   &
                            +          &
               dvr_mat(3,5) * v(5,:)   &   
                            +          &
               dvr_mat(3,6) * v(6,:)   &   
                            +          &
               dvr_mat(3,7) * v(7,:)      
  v_scr(4,:) = dvr_mat(4,1) * v(1,:)   &
                            +          &
               dvr_mat(4,2) * v(2,:)   &
                            +          &
               dvr_mat(4,3) * v(3,:)   &
                            +          &   
               dvr_mat(4,4) * v(4,:)   &
                            +          &
               dvr_mat(4,5) * v(5,:)   &   
                            +          &
               dvr_mat(4,6) * v(6,:)   &   
                            +          &
                dvr_mat(4,7) * v(7,:)      
  v_scr(5,:) = dvr_mat(5,1) * v(1,:)   &
                            +          &
               dvr_mat(5,2) * v(2,:)   &
                            +          &
               dvr_mat(5,3) * v(3,:)   &
                            +          &   
               dvr_mat(5,4) * v(4,:)   &
                            +          &
               dvr_mat(5,5) * v(5,:)   &   
                            +          &
               dvr_mat(5,6) * v(6,:)   &   
                            +          &
               dvr_mat(5,7) * v(7,:)      
  v_scr(6,:) = dvr_mat(6,1) * v(1,:)   &
                            +          &
               dvr_mat(6,2) * v(2,:)   &
                            +          &
               dvr_mat(6,3) * v(3,:)   &
                            +          &   
               dvr_mat(6,4) * v(4,:)   &
                            +          &
               dvr_mat(6,5) * v(5,:)   &   
                            +          &
               dvr_mat(6,6) * v(6,:)   &   
                            +          &
               dvr_mat(6,7) * v(7,:)      
  v_scr(7,:) = dvr_mat(7,1) * v(1,:)   &
                            +          &
               dvr_mat(7,2) * v(2,:)   &
                            +          &
               dvr_mat(7,3) * v(3,:)   &
                            +          &   
               dvr_mat(7,4) * v(4,:)   &
                            +          &
               dvr_mat(7,5) * v(5,:)   &   
                            +          &
               dvr_mat(7,6) * v(6,:)   &   
                            +          &
               dvr_mat(7,7) * v(7,:)      
!
  v(1:7,:) = v_scr(1:7,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_7_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_so_mat_v_8_d
!***begin prologue     v_so_mat_v_8_d     
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
!***end prologue       v_so_mat_v_8_d
!
  SUBROUTINE v_so_mat_v_8_d(v,             &
                            v_scr,         &
                            dvr_mat,       &
                            ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(8,8)               :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:)    = dvr_mat(1,1) * v(1,:)   &
                                 +        &
                  dvr_mat(1,2) * v(2,:)   &
                                 +        &
                  dvr_mat(1,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(1,4) * v(4,:)   &
                                 +        &
                  dvr_mat(1,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(1,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(1,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(1,8) * v(8,:)      

  v_scr(2,:) =    dvr_mat(2,1) * v(1,:)   &
                                 +        &
                  dvr_mat(2,2) * v(2,:)   &
                                 +        &
                  dvr_mat(2,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(2,4) * v(4,:)   &
                                 +        &
                  dvr_mat(2,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(2,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(2,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(2,8) * v(8,:)      

  v_scr(3,:) =    dvr_mat(3,1) * v(1,:)   &
                                 +        &
                  dvr_mat(3,2) * v(2,:)   &
                                 +        &
                  dvr_mat(3,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(3,4) * v(4,:)   &
                                 +        &
                  dvr_mat(3,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(3,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(3,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(3,8) * v(8,:)      

  v_scr(4,:) =    dvr_mat(4,1) * v(1,:)   &
                                 +        &
                  dvr_mat(4,2) * v(2,:)   &
                                 +        &
                  dvr_mat(4,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(4,4) * v(4,:)   &
                                 +        &
                  dvr_mat(4,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(4,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(4,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(4,8) * v(8,:)      

  v_scr(5,:) =    dvr_mat(5,1) * v(1,:)   &
                                 +        &
                   dvr_mat(5,2) * v(2,:)  &
                                 +        &
                  dvr_mat(5,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(5,4) * v(4,:)   &
                                 +        &
                  dvr_mat(5,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(5,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(5,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(5,8) * v(8,:)     
   
  v_scr(6,:) =    dvr_mat(6,1) * v(1,:)   &
                                 +        &
                  dvr_mat(6,2) * v(2,:)   &
                                 +        &
                  dvr_mat(6,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(6,4) * v(4,:)   &
                                 +        &
                  dvr_mat(6,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(6,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(6,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(6,8) * v(8,:)        

  v_scr(7,:) =    dvr_mat(7,1) * v(1,:)   &
                                 +        &
                  dvr_mat(7,2) * v(2,:)   &
                                 +        &
                  dvr_mat(7,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(7,4) * v(4,:)   &
                                 +        &
                  dvr_mat(7,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(7,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(7,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(7,8) * v(8,:)        

  v_scr(8,:) =    dvr_mat(8,1) * v(1,:)   &
                                 +        &
                  dvr_mat(8,2) * v(2,:)   &
                                 +        &
                  dvr_mat(8,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(8,4) * v(4,:)   &
                                 +        &
                  dvr_mat(8,5) * v(5,:)   &  
                                 +        &
                  dvr_mat(8,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(8,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(8,8) * v(8,:)        

!

  v(1:8,:) = v_scr(1:8,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_8_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_so_mat_v_9_d
!***begin prologue     v_so_mat_v_9_d     
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
!***end prologue       v_so_mat_v_9_d
!
  SUBROUTINE v_so_mat_v_9_d(v,           &
                            v_scr,       &
                            dvr_mat,     &
                            ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(9,9)               :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) =    dvr_mat(1,1) * v(1,:)   &
                                 +        &
                  dvr_mat(1,2) * v(2,:)   &
                                 +        &
                  dvr_mat(1,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(1,4) * v(4,:)   &
                                 +        &
                  dvr_mat(1,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(1,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(1,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(1,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(1,9) * v(9,:)      
  v_scr(2,:) =    dvr_mat(2,1) * v(1,:)   &
                                 +        &
                  dvr_mat(2,2) * v(2,:)   &
                                 +        &
                  dvr_mat(2,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(2,4) * v(4,:)   &
                                 +        &
                  dvr_mat(2,5) * v(5,:)   &
                                 +        &
                  dvr_mat(2,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(2,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(2,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(2,9) * v(9,:)      

  v_scr(3,:) =    dvr_mat(3,1) * v(1,:)   &
                                 +        &
                  dvr_mat(3,2) * v(2,:)   &
                                 +        &
                  dvr_mat(3,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(3,4) * v(4,:)   &
                                 +        &
                  dvr_mat(3,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(3,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(3,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(3,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(3,9) * v(9,:)      

  v_scr(4,:) =    dvr_mat(4,1) * v(1,:)   &
                                 +        &
                  dvr_mat(4,2) * v(2,:)   &
                                 +        &
                  dvr_mat(4,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(4,4) * v(4,:)   &
                                 +        &
                  dvr_mat(4,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(4,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(4,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(4,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(4,9) * v(9,:)      

  v_scr(5,:) =    dvr_mat(5,1) * v(1,:)   &
                                 +        &
                  dvr_mat(5,2) * v(2,:)   &
                                 +        &
                  dvr_mat(5,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(5,4) * v(4,:)   &
                                 +        &
                  dvr_mat(5,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(5,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(5,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(5,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(5,9) * v(9,:)      

  v_scr(6,:) =    dvr_mat(6,1) * v(1,:)   &
                                 +        &
                  dvr_mat(6,2) * v(2,:)   &
                                 +        &
                  dvr_mat(6,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(6,4) * v(4,:)   &
                                 +        &
                  dvr_mat(6,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(6,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(6,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(6,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(6,9) * v(9,:)      

  v_scr(7,:) =    dvr_mat(7,1) * v(1,:)   &
                                 +        &
                  dvr_mat(7,2) * v(2,:)   &
                                 +        &
                  dvr_mat(7,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(7,4) * v(4,:)   &
                                 +        &
                  dvr_mat(7,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(7,6) * v(6,:)   &  
                                 +        &
                  dvr_mat(7,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(7,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(7,9) * v(9,:)      

  v_scr(8,:) =    dvr_mat(8,1) * v(1,:)   &
                                 +        &
                  dvr_mat(8,2) * v(2,:)   &
                                 +        &
                  dvr_mat(8,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(8,4) * v(4,:)   &
                                 +        &
                  dvr_mat(8,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(8,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(8,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(8,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(8,9) * v(9,:)      

  v_scr(9,:) =    dvr_mat(9,1) * v(1,:)   &
                                 +        &
                  dvr_mat(9,2) * v(2,:)   &
                                 +        &
                  dvr_mat(9,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(9,4) * v(4,:)   &
                                 +        &
                  dvr_mat(9,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(9,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(9,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(9,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(9,9) * v(9,:)      
!
  v(1:9,:) = v_scr(1:9,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_9_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_so_mat_v_10_d
!***begin prologue     v_so_mat_v_10_d     
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
!***end prologue       v_so_mat_v_10_d
!
  SUBROUTINE v_so_mat_v_10_d(v,          &
                             v_scr,      &
                             dvr_mat,    &
                             ni,nj)
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  REAL*8, DIMENSION(:,:)               :: v
  REAL*8, DIMENSION(:,:)               :: v_scr
  REAL*8, DIMENSION(10,10)             :: dvr_mat
  INTEGER                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) =    dvr_mat(1,1) * v(1,:)   &
                                 +        &
                  dvr_mat(1,2) * v(2,:)   &
                                 +        &
                  dvr_mat(1,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(1,4) * v(4,:)   &
                                 +        &
                  dvr_mat(1,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(1,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(1,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(1,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(1,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(1,10) * v(10,:)    

  v_scr(2,:) =    dvr_mat(2,1) * v(1,:)   &
                                 +        &
                  dvr_mat(2,2) * v(2,:)   &
                                 +        &
                  dvr_mat(2,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(2,4) * v(4,:)   &
                                 +        &
                  dvr_mat(2,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(2,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(2,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(2,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(2,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(2,10) * v(10,:)    

  v_scr(3,:) =    dvr_mat(3,1) * v(1,:)   &
                                 +        &
                  dvr_mat(3,2) * v(2,:)   &
                                 +        &
                  dvr_mat(3,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(3,4) * v(4,:)   &
                                 +        &
                  dvr_mat(3,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(3,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(3,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(3,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(3,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(3,10) * v(10,:)    

  v_scr(4,:) =    dvr_mat(4,1) * v(1,:)   &
                                 +        &
                  dvr_mat(4,2) * v(2,:)   &
                                 +        &
                  dvr_mat(4,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(4,4) * v(4,:)   &
                                 +        &
                  dvr_mat(4,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(4,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(4,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(4,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(4,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(4,10) * v(10,:)    

  v_scr(5,:) =    dvr_mat(5,1) * v(1,:)   &
                                 +        &
                  dvr_mat(5,2) * v(2,:)  &
                                 +        &
                  dvr_mat(5,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(5,4) * v(4,:)   &
                                 +        &
                  dvr_mat(5,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(5,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(5,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(5,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(5,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(5,10) * v(10,:)    
  v_scr(6,:) =    dvr_mat(6,1) * v(1,:)   &
                                 +        &
                  dvr_mat(6,2) * v(2,:)   &
                                 +        &
                  dvr_mat(6,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(6,4) * v(4,:)   &
                                 +        &
                  dvr_mat(6,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(6,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(6,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(6,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(6,9) * v(9,:)   &  
                                 +        &
                  dvr_mat(6,10) * v(10,:)     

  v_scr(7,:) =    dvr_mat(7,1) * v(1,:)   &
                                 +        &
                  dvr_mat(7,2) * v(2,:)   &
                                 +        &
                  dvr_mat(7,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(7,4) * v(4,:)   &
                                 +        &
                  dvr_mat(7,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(7,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(7,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(7,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(7,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(7,10) * v(10,:)   

  v_scr(8,:) =    dvr_mat(8,1) * v(1,:)   &
                                 +        &
                  dvr_mat(8,2) * v(2,:)   &
                                 +        &
                  dvr_mat(8,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(8,4) * v(4,:)   &
                                 +        &
                  dvr_mat(8,5) * v(5,:)   &  
                                 +        &
                  dvr_mat(8,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(8,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(8,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(8,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(8,10) * v(10,:)      

  v_scr(9,:) =    dvr_mat(9,1) * v(1,:)   &
                                 +        &
                  dvr_mat(9,2) * v(2,:)   &
                                 +        &
                  dvr_mat(9,3) * v(3,:)   &
                                 +        &   
                  dvr_mat(9,4) * v(4,:)   &
                                 +        &
                  dvr_mat(9,5) * v(5,:)   &   
                                 +        &
                  dvr_mat(9,6) * v(6,:)   &   
                                 +        &
                  dvr_mat(9,7) * v(7,:)   &   
                                 +        &
                  dvr_mat(9,8) * v(8,:)   &   
                                 +        &
                  dvr_mat(9,9) * v(9,:)   &   
                                 +        &
                  dvr_mat(9,10) * v(10,:)    

  v_scr(10,:) =   dvr_mat(10,1) * v(1,:) &
                                   +     &
                  dvr_mat(10,2) * v(2,:) &
                                   +     &
                  dvr_mat(10,3) * v(3,:) &
                                   +     &   
                  dvr_mat(10,4) * v(4,:) &
                                   +     &
                  dvr_mat(10,5) * v(5,:) &   
                                   +     &
                  dvr_mat(10,6) * v(6,:) &   
                                   +     &
                  dvr_mat(10,7) * v(7,:) &   
                                   +     &
                  dvr_mat(10,8) * v(8,:) &   
                                   +     &
                  dvr_mat(10,9) * v(9,:) &   
                                   +     &
                  dvr_mat(10,10) * v(10,:)    
!
  v(1:10,:) = v_scr(1:10,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_10_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_gen_z
!***begin prologue     v_so_mat_v_gen_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propator matrix vector multiplies
!***                   for a general FEDVR Hamiltonian.
!***                   The parameter nj is a dummy and can take on values
!***                   consistent with any dimensional problem.  In a 1D
!***                   case nj=1, in 2D nj=nx and in 3D nj=ny*nx.
!
!***references
!***routines called    ebcxx, ambcxx, apbcxx
!***end prologue       
!
  SUBROUTINE v_so_mat_v_gen_z(v,                                   &
                              v_scr,                               &
                              dvr_mat,                             &
                              ni,nj,nk)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, nk
  INTEGER                                  :: i, k
  COMPLEX*16, DIMENSION(:,:)               :: v
  COMPLEX*16, DIMENSION(:,:)               :: v_scr
  COMPLEX*16, DIMENSION(nk,nk)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  DO i=1,nk
     DO k=1,nk
        v_scr(k,:) = v_scr(k,:) + dvr_mat(k,i) * v(i,:)     
     END DO
  END DO
!
  v(1:nk,:) = v_scr(1:nk,:)
!
END SUBROUTINE  v_so_mat_v_gen_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_2_z
!***begin prologue     v_so_mat_v_2_z     
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
!***end prologue       v_so_mat_v_2_z
!
  SUBROUTINE v_so_mat_v_2_z(v,          &
                            v_scr,      &
                            dvr_mat,    &
                            ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(:,:)               :: v
  COMPLEX*16, DIMENSION(:,:)               :: v_scr
  COMPLEX*16, DIMENSION(2,2)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)     &
                            +            &
               dvr_mat(1,2) * v(2,:)  
  v_scr(2,:) = dvr_mat(2,1) * v(1,:)     &
                            +            &
               dvr_mat(2,2) * v(2,:)     
!
  v(1:2,:) = v_scr(1:2,:)
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_2_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_so_mat_v_3_z
!***begin prologue     v_so_mat_v_3_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_3_z
!
  SUBROUTINE v_so_mat_v_3_z(v,          &
                            v_scr,      &
                            dvr_mat,    &
                            ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(:,:)               :: v
  COMPLEX*16, DIMENSION(:,:)               :: v_scr
  COMPLEX*16, DIMENSION(3,3)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)     &
                            +            &
               dvr_mat(1,2) * v(2,:)     &
                            +            &
               dvr_mat(1,3) * v(3,:)     
  v_scr(2,:) = dvr_mat(2,1) * v(1,:)     &
                            +            &
               dvr_mat(2,2) * v(2,:)     &
                            +            &
               dvr_mat(2,3) * v(3,:)     
  v_scr(3,:) = dvr_mat(3,1) * v(1,:)     &
                            +            &
               dvr_mat(3,2) * v(2,:)     &
                            +            &
               dvr_mat(3,3) * v(3,:)     
!
  v(1:3,:) = v_scr(1:3,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_3_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck v_so_mat_v_4_z
!***begin prologue     v_so_mat_v_4_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_4_z
!
  SUBROUTINE v_so_mat_v_4_z(v,          &
                            v_scr,      &
                            dvr_mat,    &
                            ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(:,:)               :: v
  COMPLEX*16, DIMENSION(:,:)               :: v_scr
  COMPLEX*16, DIMENSION(4,4)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)     &
                            +            &
               dvr_mat(1,2) * v(2,:)     &
                            +            &
               dvr_mat(1,3) * v(3,:)     &
                            +            &
               dvr_mat(1,4) * v(4,:)     

  v_scr(2,:) = dvr_mat(2,1) * v(1,:)     &
                            +            &
               dvr_mat(2,2) * v(2,:)     &
                            +            &
               dvr_mat(2,3) * v(3,:)     &
                            +            &
               dvr_mat(2,4) * v(4,:)     

  v_scr(3,:) = dvr_mat(3,1) * v(1,:)     &
                            +            &
               dvr_mat(3,2) * v(2,:)     &
                            +            &
               dvr_mat(3,3) * v(3,:)     &
                            +            &
               dvr_mat(3,4) * v(4,:)     

  v_scr(4,:) = dvr_mat(4,1) * v(1,:)     &
                            +            &
               dvr_mat(4,2) * v(2,:)     &
                            +            &
               dvr_mat(4,3) * v(3,:)     &
                            +            &
               dvr_mat(4,4) * v(4,:)     
!
  v(1:4,:) = v_scr(1:4,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_4_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_5_z
!***begin prologue     v_so_mat_v_5_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_5_z
!
  SUBROUTINE v_so_mat_v_5_z(v,          &
                            v_scr,      &
                            dvr_mat,    &
                            ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(:,:)               :: v
  COMPLEX*16, DIMENSION(:,:)               :: v_scr
  COMPLEX*16, DIMENSION(5,5)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)     &
                            +            &
               dvr_mat(1,2) * v(2,:)     &
                            +            &
               dvr_mat(1,3) * v(3,:)     &
                            +            &
               dvr_mat(1,4) * v(4,:)     &
                            +            &
               dvr_mat(1,5) * v(5,:)     

  v_scr(2,:) = dvr_mat(2,1) * v(1,:)     &
                            +            &
               dvr_mat(2,2) * v(2,:)     &
                            +            &
               dvr_mat(2,3) * v(3,:)     &
                            +            &
               dvr_mat(2,4) * v(4,:)     &
                            +            &
               dvr_mat(2,5) * v(5,:)     

  v_scr(3,:) = dvr_mat(3,1) * v(1,:)     &
                            +            &
               dvr_mat(3,2) * v(2,:)     &
                            +            &
               dvr_mat(3,3) * v(3,:)     &
                            +            &
               dvr_mat(3,4) * v(4,:)     &
                            +            &
               dvr_mat(3,5) * v(5,:)     

  v_scr(4,:) = dvr_mat(4,1) * v(1,:)     &
                            +            &
               dvr_mat(4,2) * v(2,:)     &
                            +            &
               dvr_mat(4,3) * v(3,:)     &
                            +            &
               dvr_mat(4,4) * v(4,:)     &
                            +            &
               dvr_mat(4,5) * v(5,:)     

  v_scr(5,:) = dvr_mat(5,1) * v(1,:)     &
                            +            &
               dvr_mat(5,2) * v(2,:)     &
                            +            &
               dvr_mat(5,3) * v(3,:)     &
                            +            &
               dvr_mat(5,4) * v(4,:)     &
                            +            &
               dvr_mat(5,5) * v(5,:)     

!
  v(1:5,:) = v_scr(1:5,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_5_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_6_z
!***begin prologue     v_so_mat_v_6_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_6_z
!
  SUBROUTINE v_so_mat_v_6_z(v,          &
                            v_scr,      &
                            dvr_mat,    &
                            ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(:,:)               :: v
  COMPLEX*16, DIMENSION(:,:)               :: v_scr
  COMPLEX*16, DIMENSION(6,6)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)     &
                            +            &
               dvr_mat(1,2) * v(2,:)     &
                            +            &
               dvr_mat(1,3) * v(3,:)     &
                            +            &
               dvr_mat(1,4) * v(4,:)     &
                            +            &
               dvr_mat(1,5) * v(5,:)     &
                            +            &
               dvr_mat(1,6) * v(6,:)     

  v_scr(2,:) = dvr_mat(2,1) * v(1,:)     &
                            +            &
               dvr_mat(2,2) * v(2,:)     &
                            +            &
               dvr_mat(2,3) * v(3,:)     &
                            +            &
               dvr_mat(2,4) * v(4,:)     &
                            +            &
               dvr_mat(2,5) * v(5,:)     &
                            +            &
               dvr_mat(2,6) * v(6,:)     

  v_scr(3,:) = dvr_mat(3,1) * v(1,:)     &
                            +            &
               dvr_mat(3,2) * v(2,:)     &
                            +            &
               dvr_mat(3,3) * v(3,:)     &
                            +            &
               dvr_mat(3,4) * v(4,:)     &
                            +            &
               dvr_mat(3,5) * v(5,:)     &
                            +            &
               dvr_mat(3,6) * v(6,:)     

  v_scr(4,:) = dvr_mat(4,1) * v(1,:)     &
                            +            &
               dvr_mat(4,2) * v(2,:)     &
                            +            &
               dvr_mat(4,3) * v(3,:)     &
                            +            &
               dvr_mat(4,4) * v(4,:)     &
                            +            &
               dvr_mat(4,5) * v(5,:)     &
                            +            &
               dvr_mat(4,6) * v(6,:)     

  v_scr(5,:) = dvr_mat(5,1) * v(1,:)     &
                            +            &
               dvr_mat(5,2) * v(2,:)     &
                            +            &
               dvr_mat(5,3) * v(3,:)     &
                            +            &
               dvr_mat(5,4) * v(4,:)     &
                            +            &
               dvr_mat(5,5) * v(5,:)     &
                            +            &
               dvr_mat(5,6) * v(6,:)     

  v_scr(6,:) = dvr_mat(6,1) * v(1,:)     &
                            +            &
               dvr_mat(6,2) * v(2,:)     &
                            +            &
               dvr_mat(6,3) * v(3,:)     &
                            +            &
               dvr_mat(6,4) * v(4,:)     &
                            +            &
               dvr_mat(6,5) * v(5,:)     &
                            +            &
               dvr_mat(6,6) * v(6,:)     

!
  v(1:6,:) = v_scr(1:6,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_6_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_7_z
!***begin prologue     v_so_mat_v_7_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_7_z
!
  SUBROUTINE v_so_mat_v_7_z(v,          &
                            v_scr,      &
                            dvr_mat,    &
                            ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(:,:)               :: v
  COMPLEX*16, DIMENSION(:,:)               :: v_scr
  COMPLEX*16, DIMENSION(7,7)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)     &
                            +            &
               dvr_mat(1,2) * v(2,:)     &
                            +            &
               dvr_mat(1,3) * v(3,:)     &
                            +            &
               dvr_mat(1,4) * v(4,:)     &
                            +            &
               dvr_mat(1,5) * v(5,:)     &
                            +            &
               dvr_mat(1,6) * v(6,:)     &
                            +            &
               dvr_mat(1,7) * v(7,:)     

  v_scr(2,:) = dvr_mat(2,1) * v(1,:)     &
                            +            &
               dvr_mat(2,2) * v(2,:)     &
                            +            &
               dvr_mat(2,3) * v(3,:)     &
                            +            &
               dvr_mat(2,4) * v(4,:)     &
                            +            &
               dvr_mat(2,5) * v(5,:)     &
                            +            &
               dvr_mat(2,6) * v(6,:)     &
                            +            &
               dvr_mat(2,7) * v(7,:)     

  v_scr(3,:) = dvr_mat(3,1) * v(1,:)     &
                            +            &
               dvr_mat(3,2) * v(2,:)     &
                            +            &
               dvr_mat(3,3) * v(3,:)     &
                            +            &
               dvr_mat(3,4) * v(4,:)     &
                            +            &
               dvr_mat(3,5) * v(5,:)     &
                            +            &
               dvr_mat(3,6) * v(6,:)     &
                            +            &
               dvr_mat(3,7) * v(7,:)     

  v_scr(4,:) = dvr_mat(4,1) * v(1,:)     &
                            +            &
               dvr_mat(4,2) * v(2,:)     &
                            +            &
               dvr_mat(4,3) * v(3,:)     &
                            +            &
               dvr_mat(4,4) * v(4,:)     &
                            +            &
               dvr_mat(4,5) * v(5,:)     &
                            +            &
               dvr_mat(4,6) * v(6,:)     &
                            +            &
               dvr_mat(4,7) * v(7,:)     

  v_scr(5,:) = dvr_mat(5,1) * v(1,:)     &
                            +            &
               dvr_mat(5,2) * v(2,:)     &
                            +            &
               dvr_mat(5,3) * v(3,:)     &
                            +            &
               dvr_mat(5,4) * v(4,:)     &
                            +            &
               dvr_mat(5,5) * v(5,:)     &
                            +            &
               dvr_mat(5,6) * v(6,:)     &
                            +            &
               dvr_mat(5,7) * v(7,:)     

  v_scr(6,:) = dvr_mat(6,1) * v(1,:)     &
                            +            &
               dvr_mat(6,2) * v(2,:)     &
                            +            &
               dvr_mat(6,3) * v(3,:)     &
                            +            &
               dvr_mat(6,4) * v(4,:)     &
                            +            &
               dvr_mat(6,5) * v(5,:)     &
                            +            &
               dvr_mat(6,6) * v(6,:)     &
                            +            &
               dvr_mat(6,7) * v(7,:)     

  v_scr(7,:) = dvr_mat(7,1) * v(1,:)     &
                            +            &
               dvr_mat(7,2) * v(2,:)     &
                            +            &
               dvr_mat(7,3) * v(3,:)     &
                            +            &
               dvr_mat(7,4) * v(4,:)     &
                            +            &
               dvr_mat(7,5) * v(5,:)     &
                            +            &
               dvr_mat(7,6) * v(6,:)     &
                            +            &
               dvr_mat(7,7) * v(7,:)     

!
  v(1:7,:) = v_scr(1:7,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_7_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_8_z
!***begin prologue     v_so_mat_v_8_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_8_z
!
  SUBROUTINE v_so_mat_v_8_z(v,               &
                            v_scr,           &
                            dvr_mat,         &
                            ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(:,:)               :: v
  COMPLEX*16, DIMENSION(:,:)               :: v_scr
  COMPLEX*16, DIMENSION(8,8)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)     &
                            +            &
               dvr_mat(1,2) * v(2,:)     &
                            +            &
               dvr_mat(1,3) * v(3,:)     &
                            +            &
               dvr_mat(1,4) * v(4,:)     &
                            +            &
               dvr_mat(1,5) * v(5,:)     &
                            +            &
               dvr_mat(1,6) * v(6,:)     &
                            +            &
               dvr_mat(1,7) * v(7,:)     &
                            +            &
               dvr_mat(1,8) * v(8,:)     

  v_scr(2,:) = dvr_mat(2,1) * v(1,:)     &
                            +            &
               dvr_mat(2,2) * v(2,:)     &
                            +            &
               dvr_mat(2,3) * v(3,:)     &
                            +            &
               dvr_mat(2,4) * v(4,:)     &
                            +            &
               dvr_mat(2,5) * v(5,:)     &
                            +            &
               dvr_mat(2,6) * v(6,:)     &
                            +            &
               dvr_mat(2,7) * v(7,:)     &
                            +            &
               dvr_mat(2,8) * v(8,:)     

  v_scr(3,:) = dvr_mat(3,1) * v(1,:)     &
                            +            &
               dvr_mat(3,2) * v(2,:)     &
                            +            &
               dvr_mat(3,3) * v(3,:)     &
                            +            &
               dvr_mat(3,4) * v(4,:)     &
                            +            &
               dvr_mat(3,5) * v(5,:)     &
                            +            &
               dvr_mat(3,6) * v(6,:)     &
                            +            &
               dvr_mat(3,7) * v(7,:)     &
                            +            &
               dvr_mat(3,8) * v(8,:)     

  v_scr(4,:) = dvr_mat(4,1) * v(1,:)     &
                            +            &
               dvr_mat(4,2) * v(2,:)     &
                            +            &
               dvr_mat(4,3) * v(3,:)     &
                            +            &
               dvr_mat(4,4) * v(4,:)     &
                            +            &
               dvr_mat(4,5) * v(5,:)     &
                            +            &
               dvr_mat(4,6) * v(6,:)     &
                            +            &
               dvr_mat(4,7) * v(7,:)     &
                            +            &
               dvr_mat(4,8) * v(8,:)     

  v_scr(5,:) = dvr_mat(5,1) * v(1,:)     &
                            +            &
               dvr_mat(5,2) * v(2,:)     &
                            +            &
               dvr_mat(5,3) * v(3,:)     &
                            +            &
               dvr_mat(5,4) * v(4,:)     &
                            +            &
               dvr_mat(5,5) * v(5,:)     &
                            +            &
               dvr_mat(5,6) * v(6,:)     &
                            +            &
               dvr_mat(5,7) * v(7,:)     &
                            +            &
               dvr_mat(5,8) * v(8,:)     

  v_scr(6,:) = dvr_mat(6,1) * v(1,:)     &
                            +            &
               dvr_mat(6,2) * v(2,:)     &
                            +            &
               dvr_mat(6,3) * v(3,:)     &
                            +            &
               dvr_mat(6,4) * v(4,:)     &
                            +            &
               dvr_mat(6,5) * v(5,:)     &
                            +            &
               dvr_mat(6,6) * v(6,:)     &
                            +            &
               dvr_mat(6,7) * v(7,:)     &
                            +            &
               dvr_mat(6,8) * v(8,:)     

  v_scr(7,:) = dvr_mat(7,1) * v(1,:)     &
                            +            &
               dvr_mat(7,2) * v(2,:)     &
                            +            &
               dvr_mat(7,3) * v(3,:)     &
                            +            &
               dvr_mat(7,4) * v(4,:)     &
                            +            &
               dvr_mat(7,5) * v(5,:)     &
                            +            &
               dvr_mat(7,6) * v(6,:)     &
                            +            &
               dvr_mat(7,7) * v(7,:)     &
                            +            &
               dvr_mat(7,8) * v(8,:)     

  v_scr(8,:) = dvr_mat(8,1) * v(1,:)     &
                            +            &
               dvr_mat(8,2) * v(2,:)     &
                            +            &
               dvr_mat(8,3) * v(3,:)     &
                            +            &
               dvr_mat(8,4) * v(4,:)     &
                            +            &
               dvr_mat(8,5) * v(5,:)     &
                            +            &
               dvr_mat(8,6) * v(6,:)     &
                            +            &
               dvr_mat(8,7) * v(7,:)     &
                            +            &
               dvr_mat(8,8) * v(8,:)     

!
  v(1:8,:) = v_scr(1:8,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_8_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_9_z
!***begin prologue     v_so_mat_v_9_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_9_z
!
  SUBROUTINE v_so_mat_v_9_z(v,               &
                            v_scr,           &
                            dvr_mat,         &
                            ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  COMPLEX*16, DIMENSION(:,:)               :: v
  COMPLEX*16, DIMENSION(:,:)               :: v_scr
  COMPLEX*16, DIMENSION(9,9)               :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)     &
                            +            &
               dvr_mat(1,2) * v(2,:)     &
                            +            &
               dvr_mat(1,3) * v(3,:)     &
                            +            &
               dvr_mat(1,4) * v(4,:)     &
                            +            &
               dvr_mat(1,5) * v(5,:)     &
                            +            &
               dvr_mat(1,6) * v(6,:)     &
                            +            &
               dvr_mat(1,7) * v(7,:)     &
                            +            &
               dvr_mat(1,8) * v(8,:)     &
                            +            &
               dvr_mat(1,9) * v(9,:)     

  v_scr(2,:) = dvr_mat(2,1) * v(1,:)     &
                            +            &
               dvr_mat(2,2) * v(2,:)     &
                            +            &
               dvr_mat(2,3) * v(3,:)     &
                            +            &
               dvr_mat(2,4) * v(4,:)     &
                            +            &
               dvr_mat(2,5) * v(5,:)     &
                            +            &
               dvr_mat(2,6) * v(6,:)     &
                            +            &
               dvr_mat(2,7) * v(7,:)     &
                            +            &
               dvr_mat(2,8) * v(8,:)     &
                            +            &
               dvr_mat(2,9) * v(9,:)     

  v_scr(3,:) = dvr_mat(3,1) * v(1,:)     &
                            +            &
               dvr_mat(3,2) * v(2,:)     &
                            +            &
               dvr_mat(3,3) * v(3,:)     &
                            +            &
               dvr_mat(3,4) * v(4,:)     &
                            +            &
               dvr_mat(3,5) * v(5,:)     &
                            +            &
               dvr_mat(3,6) * v(6,:)     &
                            +            &
               dvr_mat(3,7) * v(7,:)     &
                            +            &
               dvr_mat(3,8) * v(8,:)     &
                            +            &
               dvr_mat(3,9) * v(9,:)     

  v_scr(4,:) = dvr_mat(4,1) * v(1,:)     &
                            +            &
               dvr_mat(4,2) * v(2,:)     &
                            +            &
               dvr_mat(4,3) * v(3,:)     &
                            +            &
               dvr_mat(4,4) * v(4,:)     &
                            +            &
               dvr_mat(4,5) * v(5,:)     &
                            +            &
               dvr_mat(4,6) * v(6,:)     &
                            +            &
               dvr_mat(4,7) * v(7,:)     &
                            +            &
               dvr_mat(4,8) * v(8,:)     &
                            +            &
               dvr_mat(4,9) * v(9,:)     

  v_scr(5,:) = dvr_mat(5,1) * v(1,:)     &
                            +            &
               dvr_mat(5,2) * v(2,:)     &
                            +            &
               dvr_mat(5,3) * v(3,:)     &
                            +            &
               dvr_mat(5,4) * v(4,:)     &
                            +            &
               dvr_mat(5,5) * v(5,:)     &
                            +            &
               dvr_mat(5,6) * v(6,:)     &
                            +            &
               dvr_mat(5,7) * v(7,:)     &
                            +            &
               dvr_mat(5,8) * v(8,:)     &
                            +            &
               dvr_mat(5,9) * v(9,:)     

  v_scr(6,:) = dvr_mat(6,1) * v(1,:)     &
                            +            &
               dvr_mat(6,2) * v(2,:)     &
                            +            &
               dvr_mat(6,3) * v(3,:)     &
                            +            &
               dvr_mat(6,4) * v(4,:)     &
                            +            &
               dvr_mat(6,5) * v(5,:)     &
                            +            &
               dvr_mat(6,6) * v(6,:)     &
                            +            &
               dvr_mat(6,7) * v(7,:)     &
                            +            &
               dvr_mat(6,8) * v(8,:)     &
                            +            &
               dvr_mat(6,9) * v(9,:)     

  v_scr(7,:) = dvr_mat(7,1) * v(1,:)     &
                            +            &
               dvr_mat(7,2) * v(2,:)     &
                            +            &
               dvr_mat(7,3) * v(3,:)     &
                            +            &
               dvr_mat(7,4) * v(4,:)     &
                            +            &
               dvr_mat(7,5) * v(5,:)     &
                            +            &
               dvr_mat(7,6) * v(6,:)     &
                            +            &
               dvr_mat(7,7) * v(7,:)     &
                            +            &
               dvr_mat(7,8) * v(8,:)     &
                            +            &
               dvr_mat(7,9) * v(9,:)     

  v_scr(8,:) = dvr_mat(8,1) * v(1,:)     &
                            +            &
               dvr_mat(8,2) * v(2,:)     &
                            +            &
               dvr_mat(8,3) * v(3,:)     &
                            +            &
               dvr_mat(8,4) * v(4,:)     &
                            +            &
               dvr_mat(8,5) * v(5,:)     &
                            +            &
               dvr_mat(8,6) * v(6,:)     &
                            +            &
               dvr_mat(8,7) * v(7,:)     &
                            +            &
               dvr_mat(8,8) * v(8,:)     &
                            +            &
               dvr_mat(8,9) * v(9,:)     

  v_scr(9,:) = dvr_mat(9,1) * v(1,:)     &
                            +            &
               dvr_mat(9,2) * v(2,:)     &
                            +            &
               dvr_mat(9,3) * v(3,:)     &
                            +            &
               dvr_mat(9,4) * v(4,:)     &
                            +            &
               dvr_mat(9,5) * v(5,:)     &
                            +            &
               dvr_mat(9,6) * v(6,:)     &
                            +            &
               dvr_mat(9,7) * v(7,:)     &
                            +            &
               dvr_mat(9,8) * v(8,:)     &
                            +            &
               dvr_mat(9,9) * v(9,:)     
!
  v(1:9,:) = v_scr(1:9,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_9_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_so_mat_v_10_z
!***begin prologue     v_so_mat_v_10_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multiplies
!***                   for the special case of 3*3 matrices.
!
!***references
!***routines called
!***end prologue       v_so_mat_v_10_z
!
  SUBROUTINE v_so_mat_v_10_z(v,               &
                             v_scr,           &
                             dvr_mat,         &
                             ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                              :: ni, nj
  COMPLEX*16, DIMENSION(:,:)               :: v
  COMPLEX*16, DIMENSION(:,:)               :: v_scr
  COMPLEX*16, DIMENSION(10,10)             :: dvr_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  v_scr(1,:) = dvr_mat(1,1) * v(1,:)     &
                            +            &
               dvr_mat(1,2) * v(2,:)     &
                            +            &
               dvr_mat(1,3) * v(3,:)     &
                            +            &
               dvr_mat(1,4) * v(4,:)     &
                            +            &
               dvr_mat(1,5) * v(5,:)     &
                            +            &
               dvr_mat(1,6) * v(6,:)     &
                            +            &
               dvr_mat(1,7) * v(7,:)     &
                            +            &
               dvr_mat(1,8) * v(8,:)     &
                            +            &
               dvr_mat(1,9) * v(9,:)     &
                            +            &
               dvr_mat(1,10) * v(10,:)   

  v_scr(2,:) = dvr_mat(2,1) * v(1,:)     &
                            +            &
               dvr_mat(2,2) * v(2,:)     &
                            +            &
               dvr_mat(2,3) * v(3,:)     &
                            +            &
               dvr_mat(2,4) * v(4,:)     &
                            +            &
               dvr_mat(2,5) * v(5,:)     &
                            +            &
               dvr_mat(2,6) * v(6,:)     &
                            +            &
               dvr_mat(2,7) * v(7,:)     &
                            +            &
               dvr_mat(2,8) * v(8,:)     &
                            +            &
               dvr_mat(2,9) * v(9,:)     &
                            +            &
               dvr_mat(2,10) * v(10,:)

  v_scr(3,:) = dvr_mat(3,1) * v(1,:)     &
                            +            &
               dvr_mat(3,2) * v(2,:)     &
                            +            &
               dvr_mat(3,3) * v(3,:)     &
                            +            &
               dvr_mat(3,4) * v(4,:)     &
                            +            &
               dvr_mat(3,5) * v(5,:)     &
                            +            &
               dvr_mat(3,6) * v(6,:)     &
                            +            &
               dvr_mat(3,7) * v(7,:)     &
                            +            &
               dvr_mat(3,8) * v(8,:)     &
                            +            &
               dvr_mat(3,9) * v(9,:)     &
                            +            &
               dvr_mat(3,10) * v(10,:)   

  v_scr(4,:) = dvr_mat(4,1) * v(1,:)     &
                            +            &
               dvr_mat(4,2) * v(2,:)     &
                            +            &
               dvr_mat(4,3) * v(3,:)     &
                            +            &
               dvr_mat(4,4) * v(4,:)     &
                            +            &
               dvr_mat(4,5) * v(5,:)     &
                            +            &
               dvr_mat(4,6) * v(6,:)     &
                            +            &
               dvr_mat(4,7) * v(7,:)     &
                            +            &
               dvr_mat(4,8) * v(8,:)     &
                            +            &
               dvr_mat(4,9) * v(9,:)     &
                            +            &
               dvr_mat(4,10) * v(10,:)   

  v_scr(5,:) = dvr_mat(5,1) * v(1,:)     &
                            +            &
               dvr_mat(5,2) * v(2,:)     &
                            +            &
               dvr_mat(5,3) * v(3,:)     &
                            +            &
               dvr_mat(5,4) * v(4,:)     &
                            +            &
               dvr_mat(5,5) * v(5,:)     &
                            +            &
               dvr_mat(5,6) * v(6,:)     &
                            +            &
               dvr_mat(5,7) * v(7,:)     &
                            +            &
               dvr_mat(5,8) * v(8,:)     &
                            +            &
               dvr_mat(5,9) * v(9,:)     &
                            +            &
               dvr_mat(5,10) * v(10,:)   

  v_scr(6,:) = dvr_mat(6,1) * v(1,:)     &
                            +            &
               dvr_mat(6,2) * v(2,:)     &
                            +            &
               dvr_mat(6,3) * v(3,:)     &
                            +            &
               dvr_mat(6,4) * v(4,:)     &
                            +            &
               dvr_mat(6,5) * v(5,:)     &
                            +            &
               dvr_mat(6,6) * v(6,:)     &
                            +            &
               dvr_mat(6,7) * v(7,:)     &
                            +            &
               dvr_mat(6,8) * v(8,:)     &
                            +            &
               dvr_mat(6,9) * v(9,:)     &
                            +            &
               dvr_mat(6,10) * v(10,:)   

  v_scr(7,:) = dvr_mat(7,1) * v(1,:)     &
                            +            &
               dvr_mat(7,2) * v(2,:)     &
                            +            &
               dvr_mat(7,3) * v(3,:)     &
                            +            &
               dvr_mat(7,4) * v(4,:)     &
                            +            &
               dvr_mat(7,5) * v(5,:)     &
                            +            &
               dvr_mat(7,6) * v(6,:)     &
                            +            &
               dvr_mat(7,7) * v(7,:)     &
                            +            &
               dvr_mat(7,8) * v(8,:)     &
                            +            &
               dvr_mat(7,9) * v(9,:)     &
                            +            &
               dvr_mat(7,10) * v(10,:)   

  v_scr(8,:) = dvr_mat(8,1) * v(1,:)     &
                            +            &
               dvr_mat(8,2) * v(2,:)     &
                            +            &
               dvr_mat(8,3) * v(3,:)     &
                            +            &
               dvr_mat(8,4) * v(4,:)     &
                            +            &
               dvr_mat(8,5) * v(5,:)     &
                            +            &
               dvr_mat(8,6) * v(6,:)     &
                            +            &
               dvr_mat(8,7) * v(7,:)     &
                            +            &
               dvr_mat(8,8) * v(8,:)     &
                            +            &
               dvr_mat(8,9) * v(9,:)     &
                            +            &
               dvr_mat(8,10) * v(10,:)   

  v_scr(9,:) = dvr_mat(9,1) * v(1,:)     &
                            +            &
               dvr_mat(9,2) * v(2,:)     &
                            +            &
               dvr_mat(9,3) * v(3,:)     &
                            +            &
               dvr_mat(9,4) * v(4,:)     &
                            +            &
               dvr_mat(9,5) * v(5,:)     &
                            +            &
               dvr_mat(9,6) * v(6,:)     &
                            +            &
               dvr_mat(9,7) * v(7,:)     &
                            +            &
               dvr_mat(9,8) * v(8,:)     &
                            +            &
               dvr_mat(9,9) * v(9,:)     &
                            +            &
               dvr_mat(9,10) * v(10,:)   

  v_scr(10,:) = dvr_mat(10,1) * v(1,:)   &
                              +          &
                dvr_mat(10,2) * v(2,:)   &
                              +          &
                dvr_mat(10,3) * v(3,:)   &
                              +          &
                dvr_mat(10,4) * v(4,:)   &
                              +          &
                dvr_mat(10,5) * v(5,:)   &
                              +          &
                dvr_mat(10,6) * v(6,:)   &
                              +          &
                dvr_mat(10,7) * v(7,:)   &
                              +          &
                dvr_mat(10,8) * v(8,:)   &
                              +          &
                dvr_mat(10,9) * v(9,:)   &
                              +          &
                dvr_mat(10,10) * v(10,:)   
!
  v(1:10,:) = v_scr(1:10,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_so_mat_v_10_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               END       MODULE v_out_so_mat_v_in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

