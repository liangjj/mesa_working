!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                 MODULE v_m_v
!**begin prologue     v_m_v
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       v_m_v
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                 CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_m_v_gen
!***begin prologue     v_m_v_gen     
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
  SUBROUTINE v_m_v_gen(v,          &
                       v_scr,      &
                       exp_t_mat,  &
                       ni,nj,nk)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, nk
  INTEGER                                  :: i, j, k
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(nk,nk)                 :: exp_t_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   v_scr(1:nk,1:nj) = 0.d0
   DO i=1,nk
      DO k=1,nk
         DO j=1,nj
            v_scr(k,j) = v_scr(k,j) + exp_t_mat(k,i) * v(i,j) 
         END DO
      END DO
   END DO
!
!       Copy the temporary vector back to the input vector.
!
  v(1:nk,:) = v_scr(1:nk,:)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE v_m_v_gen
!
!
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
  SUBROUTINE v_m_v_2(v,                   &
                     v_scr,               &
                     exp_t_mat,           &
                     ni,nj)
  USE io
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(2,2)                   :: exp_t_mat
  INTEGER                                  :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  DO i=1,nj
     v_scr(1,i) = exp_t_mat(1,1) * v(1,i) + exp_t_mat(1,2) * v(2,i)
     v_scr(2,i) = exp_t_mat(2,1) * v(1,i) + exp_t_mat(2,2) * v(2,i)
  END DO
  v(1:2,:) = v_scr(1:2,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_m_v_2
!
!
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
  SUBROUTINE v_m_v_3(v,              &
                     v_scr,          &
                     exp_t_mat,      &
                     ni,nj)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(3,3)                   :: exp_t_mat
  INTEGER                                  :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  DO i=1,nj
     v_scr(1,i) = exp_t_mat(1,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(1,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(1,3) * v(3,i)   
     v_scr(2,i) = exp_t_mat(2,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(2,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(2,3) * v(3,i)   
     v_scr(3,i) = exp_t_mat(3,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(3,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(3,3) * v(3,i)   
  END DO
  v(1:3,:) = v_scr(1:3,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_m_v_3
!
!
!deck v_m_v_4
!***begin prologue     v_m_v_4     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multipleis
!***                   for the special case of 4*4 matrices.
!
!***references
!***routines called
!***end prologue       v_m_v_4
!
  SUBROUTINE v_m_v_4(v,             &
                     v_scr,         &
                     exp_t_mat,     &
                     ni,nj)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(4,4)                   :: exp_t_mat
  INTEGER                                  :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  DO i=1,nj
     v_scr(1,i) = exp_t_mat(1,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(1,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(1,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(1,4) * v(4,i)   

     v_scr(2,i) = exp_t_mat(2,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(2,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(2,3) * v(3,i)   &
                                 +          &
                  exp_t_mat(2,4) * v(4,i)   

     v_scr(3,i) = exp_t_mat(3,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(3,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(3,3) * v(3,i)   &
                                 +          &
                  exp_t_mat(3,4) * v(4,i)   

     v_scr(4,i) = exp_t_mat(4,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(4,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(4,3) * v(3,i)   &
                                 +          &
                  exp_t_mat(4,4) * v(4,i)   
  END DO
  v(1:4,:) = v_scr(1:4,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_m_v_4
!
!
!deck v_m_v_5
!***begin prologue     v_m_v_5     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        Compute the propagator matrix vector multipleis
!***                   for the special case of 5*5 matrices.
!
!***references
!***routines called
!***end prologue       v_m_v_5
!
  SUBROUTINE v_m_v_5(v,            &
                     v_scr,        &
                     exp_t_mat,    &
                     ni,nj)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(5,5)                   :: exp_t_mat
  INTEGER                                  :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  DO i=1,nj
     v_scr(1,i) = exp_t_mat(1,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(1,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(1,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(1,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(1,5) * v(5,i)      
     v_scr(2,i) = exp_t_mat(2,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(2,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(2,3) * v(3,i)   &
                                 +          &
                  exp_t_mat(2,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(2,5) * v(5,i)   
     v_scr(3,i) = exp_t_mat(3,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(3,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(3,3) * v(3,i)   &
                                 +          &
                  exp_t_mat(3,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(3,5) * v(5,i)   
     v_scr(4,i) = exp_t_mat(4,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(4,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(4,3) * v(3,i)   &
                                 +          &
                  exp_t_mat(4,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(4,5) * v(5,i)    
     v_scr(5,i) = exp_t_mat(5,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(5,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(5,3) * v(3,i)   &
                                 +          &
                  exp_t_mat(5,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(5,5) * v(5,i)    
  END DO
  v(1:5,:) = v_scr(1:5,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_m_v_5
!
!
!deck v_m_v_6
!***begin prologue     v_m_v_6     
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
!***end prologue       v_m_v_6
!
  SUBROUTINE v_m_v_6(v,           &    
                     v_scr,       &
                     exp_t_mat,   &
                     ni,nj)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(6,6)                   :: exp_t_mat
  INTEGER                                  :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  DO i=1,nj
     v_scr(1,i) = exp_t_mat(1,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(1,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(1,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(1,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(1,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(1,6) * v(6,i)      

     v_scr(2,i) = exp_t_mat(2,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(2,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(2,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(2,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(2,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(2,6) * v(6,i)      
     v_scr(3,i) = exp_t_mat(3,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(3,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(3,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(3,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(3,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(3,6) * v(6,i)      
     v_scr(4,i) = exp_t_mat(4,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(4,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(4,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(4,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(4,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(4,6) * v(6,i)      
     v_scr(5,i) = exp_t_mat(5,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(5,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(5,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(5,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(5,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(5,6) * v(6,i)      

     v_scr(6,i) = exp_t_mat(6,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(6,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(6,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(6,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(6,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(6,6) * v(6,i)      
  END DO
  v(1:6,:) = v_scr(1:6,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_m_v_6
!
!
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
  SUBROUTINE v_m_v_7(v,           &    
                     v_scr,       &
                     exp_t_mat,   &
                     ni,nj)
                     
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(7,7)                   :: exp_t_mat
  INTEGER                                  :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  DO i=1,nj
     v_scr(1,i) = exp_t_mat(1,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(1,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(1,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(1,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(1,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(1,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(1,7) * v(7,i)      

     v_scr(2,i) = exp_t_mat(2,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(2,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(2,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(2,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(2,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(2,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(2,7) * v(7,i)      
     v_scr(3,i) = exp_t_mat(3,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(3,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(3,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(3,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(3,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(3,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(3,7) * v(7,i)      
     v_scr(4,i) = exp_t_mat(4,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(4,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(4,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(4,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(4,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(4,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(4,7) * v(7,i)      
     v_scr(5,i) = exp_t_mat(5,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(5,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(5,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(5,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(5,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(5,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(5,7) * v(7,i)      
     v_scr(6,i) = exp_t_mat(6,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(6,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(6,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(6,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(6,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(6,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(6,7) * v(7,i)      
     v_scr(7,i) = exp_t_mat(7,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(7,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(7,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(7,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(7,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(7,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(7,7) * v(7,i)      
  END DO
  v(1:7,:) = v_scr(1:7,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_m_v_7
!
!
!deck v_m_v_8
!***begin prologue     v_m_v_8     
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
!***end prologue       v_m_v_8
!
  SUBROUTINE v_m_v_8(v,              &
                     v_scr,          &
                     exp_t_mat,      &
                     ni,nj)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(8,8)                   :: exp_t_mat
  INTEGER                                  :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  DO i=1,nj
     v_scr(1,i) = exp_t_mat(1,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(1,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(1,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(1,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(1,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(1,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(1,7) * v(7,i)   &   
                                 +          &
                  exp_t_mat(1,8) * v(8,i)      

     v_scr(2,i) = exp_t_mat(2,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(2,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(2,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(2,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(2,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(2,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(2,7) * v(7,i)   &   
                                 +          &
                  exp_t_mat(2,8) * v(8,i)      

     v_scr(3,i) = exp_t_mat(3,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(3,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(3,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(3,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(3,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(3,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(3,7) * v(7,i)   &   
                                 +          &
                  exp_t_mat(3,8) * v(8,i)      

     v_scr(4,i) = exp_t_mat(4,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(4,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(4,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(4,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(4,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(4,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(4,7) * v(7,i)   &   
                                 +          &
                  exp_t_mat(4,8) * v(8,i)      

     v_scr(5,i) = exp_t_mat(5,1) * v(1,i)     &
                                 +            &
                   exp_t_mat(5,2) * v(2,i)    &
                                 +            &
                  exp_t_mat(5,3) * v(3,i)     &
                                 +            &   
                  exp_t_mat(5,4) * v(4,i)     &
                                 +            &
                  exp_t_mat(5,5) * v(5,i)     &   
                                 +            &
                  exp_t_mat(5,6) * v(6,i)     &   
                                 +            &
                  exp_t_mat(5,7) * v(7,i)     &   
                                 +            &
                  exp_t_mat(5,8) * v(8,i)     
   
     v_scr(6,i) = exp_t_mat(6,1) * v(1,i)     &
                                 +            &
                  exp_t_mat(6,2) * v(2,i)     &
                                 +            &
                  exp_t_mat(6,3) * v(3,i)     &
                                 +            &   
                  exp_t_mat(6,4) * v(4,i)     &
                                 +            &
                  exp_t_mat(6,5) * v(5,i)     &   
                                 +            &
                  exp_t_mat(6,6) * v(6,i)     &   
                                 +            &
                  exp_t_mat(6,7) * v(7,i)     &   
                                 +            &
                  exp_t_mat(6,8) * v(8,i)        

     v_scr(7,i) = exp_t_mat(7,1) * v(1,i)     &
                                 +            &
                  exp_t_mat(7,2) * v(2,i)     &
                                 +            &
                  exp_t_mat(7,3) * v(3,i)     &
                                 +            &   
                  exp_t_mat(7,4) * v(4,i)     &
                                 +            &
                  exp_t_mat(7,5) * v(5,i)     &   
                                 +            &
                  exp_t_mat(7,6) * v(6,i)     &   
                                 +            &
                  exp_t_mat(7,7) * v(7,i)     &   
                                 +            &
                  exp_t_mat(7,8) * v(8,i)        

     v_scr(8,i) = exp_t_mat(8,1) * v(1,i)     &
                                 +            &
                  exp_t_mat(8,2) * v(2,i)     &
                                 +            &
                  exp_t_mat(8,3) * v(3,i)     &
                                 +            &   
                  exp_t_mat(8,4) * v(4,i)     &
                                 +            &
                  exp_t_mat(8,5) * v(5,i)     &  
                                 +            &
                  exp_t_mat(8,6) * v(6,i)     &   
                                 +            &
                  exp_t_mat(8,7) * v(7,i)     &   
                                 +            &
                  exp_t_mat(8,8) * v(8,i)        

  END DO
  v(1:8,:) = v_scr(1:8,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_m_v_8
!
!
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
  SUBROUTINE v_m_v_9(v,             &
                     v_scr,         &
                     exp_t_mat,     &
                     ni,nj)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(9,9)                   :: exp_t_mat
  INTEGER                                  :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  DO i=1,nj
     v_scr(1,i) = exp_t_mat(1,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(1,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(1,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(1,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(1,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(1,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(1,7) * v(7,i)   &   
                                 +          &
                  exp_t_mat(1,8) * v(8,i)   &   
                                 +          &
                  exp_t_mat(1,9) * v(9,i)      
     v_scr(2,i) = exp_t_mat(2,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(2,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(2,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(2,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(2,5) * v(5,i)   &
                                 +          &
                  exp_t_mat(2,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(2,7) * v(7,i)   &   
                                 +          &
                  exp_t_mat(2,8) * v(8,i)   &   
                                 +          &
                  exp_t_mat(2,9) * v(9,i)      

     v_scr(3,i) = exp_t_mat(3,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(3,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(3,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(3,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(3,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(3,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(3,7) * v(7,i)   &   
                                 +          &
                  exp_t_mat(3,8) * v(8,i)   &   
                                 +          &
                  exp_t_mat(3,9) * v(9,i)      

     v_scr(4,i) = exp_t_mat(4,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(4,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(4,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(4,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(4,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(4,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(4,7) * v(7,i)   &   
                                 +          &
                  exp_t_mat(4,8) * v(8,i)   &   
                                 +          &
                  exp_t_mat(4,9) * v(9,i)      

     v_scr(5,i) = exp_t_mat(5,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(5,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(5,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(5,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(5,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(5,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(5,7) * v(7,i)   &   
                                 +          &
                  exp_t_mat(5,8) * v(8,i)   &   
                                 +          &
                  exp_t_mat(5,9) * v(9,i)      

     v_scr(6,i) = exp_t_mat(6,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(6,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(6,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(6,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(6,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(6,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(6,7) * v(7,i)   &   
                                 +          &
                  exp_t_mat(6,8) * v(8,i)   &   
                                 +          &
                  exp_t_mat(6,9) * v(9,i)      

     v_scr(7,i) = exp_t_mat(7,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(7,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(7,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(7,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(7,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(7,6) * v(6,i)   &  
                                 +          &
                  exp_t_mat(7,7) * v(7,i)   &   
                                 +          &
                  exp_t_mat(7,8) * v(8,i)   &   
                                 +          &
                  exp_t_mat(7,9) * v(9,i)      

     v_scr(8,i) = exp_t_mat(8,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(8,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(8,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(8,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(8,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(8,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(8,7) * v(7,i)   &   
                                 +          &
                  exp_t_mat(8,8) * v(8,i)   &   
                                 +          &
                  exp_t_mat(8,9) * v(9,i)      

     v_scr(9,i) = exp_t_mat(9,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(9,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(9,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(9,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(9,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(9,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(9,7) * v(7,i)   &   
                                 +          &
                  exp_t_mat(9,8) * v(8,i)   &   
                                 +          &
                  exp_t_mat(9,9) * v(9,i)      
  END DO
  v(1:9,:) = v_scr(1:9,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_m_v_9
!
!
!
!
!deck v_m_v_10
!***begin prologue     v_m_v_10     
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
!***end prologue       v_m_v_10
!
  SUBROUTINE v_m_v_10(v,           &
                      v_scr,       &
                      exp_t_mat,   &
                      ni,nj)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  REAL*8, DIMENSION(10,10)                 :: exp_t_mat
  INTEGER                                  :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
  DO i=1,nj
     v_scr(1,i) = exp_t_mat(1,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(1,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(1,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(1,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(1,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(1,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(1,7) * v(7,i)   &   
                                 +          &
                  exp_t_mat(1,8) * v(8,i)   &   
                                 +          &
                  exp_t_mat(1,9) * v(9,i)   &   
                                 +          &
                  exp_t_mat(1,10) * v(10,i)    

     v_scr(2,i) = exp_t_mat(2,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(2,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(2,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(2,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(2,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(2,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(2,7) * v(7,i)   &   
                                 +          &
                  exp_t_mat(2,8) * v(8,i)   &   
                                 +          &
                  exp_t_mat(2,9) * v(9,i)   &   
                                 +          &
                  exp_t_mat(2,10) * v(10,i)    

     v_scr(3,i) = exp_t_mat(3,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(3,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(3,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(3,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(3,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(3,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(3,7) * v(7,i)   &   
                                 +          &
                  exp_t_mat(3,8) * v(8,i)   &   
                                 +          &
                  exp_t_mat(3,9) * v(9,i)   &   
                                 +          &
                  exp_t_mat(3,10) * v(10,i)    

     v_scr(4,i) = exp_t_mat(4,1) * v(1,i)   &
                                 +          &
                  exp_t_mat(4,2) * v(2,i)   &
                                 +          &
                  exp_t_mat(4,3) * v(3,i)   &
                                 +          &   
                  exp_t_mat(4,4) * v(4,i)   &
                                 +          &
                  exp_t_mat(4,5) * v(5,i)   &   
                                 +          &
                  exp_t_mat(4,6) * v(6,i)   &   
                                 +          &
                  exp_t_mat(4,7) * v(7,i)   &   
                                 +          &
                  exp_t_mat(4,8) * v(8,i)   &   
                                 +          &
                  exp_t_mat(4,9) * v(9,i)   &   
                                 +          &
                  exp_t_mat(4,10) * v(10,i)    

     v_scr(5,i) = exp_t_mat(5,1) * v(1,i)     &
                                 +            &
                   exp_t_mat(5,2) * v(2,i)    &
                                 +            &
                  exp_t_mat(5,3) * v(3,i)     &
                                 +            &   
                  exp_t_mat(5,4) * v(4,i)     &
                                 +            &
                  exp_t_mat(5,5) * v(5,i)     &   
                                 +            &
                  exp_t_mat(5,6) * v(6,i)     &   
                                 +            &
                  exp_t_mat(5,7) * v(7,i)     &   
                                 +            &
                  exp_t_mat(5,8) * v(8,i)     &   
                                 +            &
                  exp_t_mat(5,9) * v(9,i)     &   
                                 +            &
                  exp_t_mat(5,10) * v(10,i)    
     v_scr(6,i) = exp_t_mat(6,1) * v(1,i)     &
                                 +            &
                  exp_t_mat(6,2) * v(2,i)     &
                                 +            &
                  exp_t_mat(6,3) * v(3,i)     &
                                 +            &   
                  exp_t_mat(6,4) * v(4,i)     &
                                 +            &
                  exp_t_mat(6,5) * v(5,i)     &   
                                 +            &
                  exp_t_mat(6,6) * v(6,i)     &   
                                 +            &
                  exp_t_mat(6,7) * v(7,i)     &   
                                 +            &
                  exp_t_mat(6,8) * v(8,i)     &   
                                 +            &
                  exp_t_mat(6,9) * v(9,i)     &  
                                 +            &
                  exp_t_mat(6,10) * v(10,i)     

     v_scr(7,i) = exp_t_mat(7,1) * v(1,i)     &
                                 +            &
                  exp_t_mat(7,2) * v(2,i)     &
                                 +            &
                  exp_t_mat(7,3) * v(3,i)     &
                                 +            &   
                  exp_t_mat(7,4) * v(4,i)     &
                                 +            &
                  exp_t_mat(7,5) * v(5,i)     &   
                                 +            &
                  exp_t_mat(7,6) * v(6,i)     &   
                                 +            &
                  exp_t_mat(7,7) * v(7,i)     &   
                                 +            &
                  exp_t_mat(7,8) * v(8,i)     &   
                                 +            &
                  exp_t_mat(7,9) * v(9,i)     &   
                                 +            &
                  exp_t_mat(7,10) * v(10,i)    

     v_scr(8,i) = exp_t_mat(8,1) * v(1,i)     &
                                 +            &
                  exp_t_mat(8,2) * v(2,i)     &
                                 +            &
                  exp_t_mat(8,3) * v(3,i)     &
                                 +            &   
                  exp_t_mat(8,4) * v(4,i)     &
                                 +            &
                  exp_t_mat(8,5) * v(5,i)     &  
                                 +            &
                  exp_t_mat(8,6) * v(6,i)     &   
                                 +            &
                  exp_t_mat(8,7) * v(7,i)     &   
                                 +            &
                  exp_t_mat(8,8) * v(8,i)     &   
                                 +            &
                  exp_t_mat(8,9) * v(9,i)     &   
                                 +            &
                  exp_t_mat(8,10) * v(10,i)      

     v_scr(9,i) = exp_t_mat(9,1) * v(1,i)     &
                                 +            &
                  exp_t_mat(9,2) * v(2,i)     &
                                 +            &
                  exp_t_mat(9,3) * v(3,i)     &
                                 +            &   
                  exp_t_mat(9,4) * v(4,i)     &
                                 +            &
                  exp_t_mat(9,5) * v(5,i)     &   
                                 +            &
                  exp_t_mat(9,6) * v(6,i)     &   
                                 +            &
                  exp_t_mat(9,7) * v(7,i)     &   
                                 +            &
                  exp_t_mat(9,8) * v(8,i)     &   
                                 +            &
                  exp_t_mat(9,9) * v(9,i)     &   
                                 +            &
                  exp_t_mat(9,10) * v(10,i)    

     v_scr(10,i) = exp_t_mat(10,1) * v(1,i)   &
                                   +          &
                   exp_t_mat(10,2) * v(2,i)   &
                                   +          &
                   exp_t_mat(10,3) * v(3,i)   &
                                   +          &   
                   exp_t_mat(10,4) * v(4,i)   &
                                   +          &
                   exp_t_mat(10,5) * v(5,i)   &   
                                   +          &
                   exp_t_mat(10,6) * v(6,i)   &   
                                   +          &
                   exp_t_mat(10,7) * v(7,i)   &   
                                   +          &
                   exp_t_mat(10,8) * v(8,i)   &   
                                   +          &
                   exp_t_mat(10,9) * v(9,i)   &   
                                   +          &
                   exp_t_mat(10,10) * v(10,i)    
  END DO
  v(1:10,:) = v_scr(1:10,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE v_m_v_10
END MODULE v_m_v
