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

  CONTAINS
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
  SUBROUTINE v_m_v_gen(real_v,imag_v,              &
                       real_v_scr,imag_v_scr,      &
                       cosine_t_mat,sine_t_mat,    &
                       ni,nj,nk)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, nk
  INTEGER                                  :: i, j, k
  REAL*8, DIMENSION(ni,nj)                 :: real_v, imag_v
  REAL*8, DIMENSION(ni,nj)                 :: real_v_scr, imag_v_scr
  REAL*8, DIMENSION(nk,nk)                 :: cosine_t_mat, sine_t_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                       General Code
!
!                 Real(V_scr) = cos_t_mat * Real(V) 
!
!  CALL ebcxx(real_v_scr,cosine_t_mat,real_v,nk,nk,nj,ni,nk,ni)
!
!                 Real(V_scr) = Real(V_scr) - sin_t_mat * Imag(V)
!
!  CALL ambcxx(real_v_scr,sine_t_mat,imag_v,nk,nk,nj,ni,nk,ni)
!
!                 Imag(V_scr) = cos_t_mat * Imag(V) 
!
!  CALL ebcxx(imag_v_scr,cosine_t_mat,imag_v,nk,nk,nj,ni,nk,ni)
!
!                 Imag(V_scr) = Imag(V_scr) + sine_t_mat * Real(V) 
!
!  CALL apbcxx(imag_v_scr,sine_t_mat,real_v,nk,nk,nj,ni,nk,ni)
   real_v_scr(1:nk,1:nj) = 0.d0
   imag_v_scr(1:nk,1:nj) = 0.d0
   DO i=1,nk
      DO k=1,nk
         DO j=1,nj
            real_v_scr(k,j) = real_v_scr(k,j) + cosine_t_mat(k,i) * real_v(i,j) &
                                              -   sine_t_mat(k,i) * imag_v(i,j)
            imag_v_scr(k,j) = imag_v_scr(k,j) + cosine_t_mat(k,i) * imag_v(i,j) &
                                              +   sine_t_mat(k,i) * real_v(i,j)
         END DO
      END DO
   END DO
!
!       Copy the temporary vector back to the input vector.
!
  real_v(1:nk,:) = real_v_scr(1:nk,:)
  imag_v(1:nk,:) = imag_v_scr(1:nk,:)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE v_m_v_gen
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
  SUBROUTINE v_m_v_10(real_v,imag_v,real_v_scr,imag_v_scr, &
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
                                         +               &
                       cosine_t_mat(1,10) * real_v(10,i) &   
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
                       sine_t_mat(1,9)   * imag_v(9,i)   &
                                         -               &
                       sine_t_mat(1,10)  * imag_v(10,i)   
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
                                         +               &
                       cosine_t_mat(2,10) * real_v(10,i)   &   
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
                       sine_t_mat(2,9)   * imag_v(9,i)   &
                                         -               &
                       sine_t_mat(2,10)   * imag_v(10,i)   
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
END SUBROUTINE v_m_v_10
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
  SUBROUTINE v_m_v_4(real_v,imag_v,real_v_scr,imag_v_scr, &
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
  DO i=1,nj
     real_v_scr(1,i) = cosine_t_mat(1,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(1,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(1,3) * real_v(3,i)   &
                                         +               &   
                       cosine_t_mat(1,4) * real_v(4,i)   &
                                         -               &
                       sine_t_mat(1,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(1,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(1,3)   * imag_v(3,i)   &
                                         -               &
                       sine_t_mat(1,4)   * imag_v(4,i)  
     real_v_scr(2,i) = cosine_t_mat(2,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(2,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(2,3) * real_v(3,i)   &
                                         +               &
                       cosine_t_mat(2,4) * real_v(4,i)   &
                                         -               &
                       sine_t_mat(2,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(2,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(2,3)   * imag_v(3,i)   &
                                         -               &
                       sine_t_mat(2,4)   * imag_v(4,i)  
     real_v_scr(3,i) = cosine_t_mat(3,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(3,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(3,3) * real_v(3,i)   &
                                         +               &
                       cosine_t_mat(3,4) * real_v(4,i)   &
                                         -               &
                       sine_t_mat(3,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(3,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(3,3)   * imag_v(3,i)   &
                                         -               &
                       sine_t_mat(3,4)   * imag_v(4,i)  
     real_v_scr(4,i) = cosine_t_mat(4,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(4,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(4,3) * real_v(3,i)   &
                                         +               &
                       cosine_t_mat(4,4) * real_v(4,i)   &
                                         -               &
                       sine_t_mat(4,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(4,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(4,3)   * imag_v(3,i)   &
                                         -               &
                       sine_t_mat(4,4)   * imag_v(4,i)  
     imag_v_scr(1,i) = cosine_t_mat(1,1) * imag_v(1,i)   &
                                         +               &
                       cosine_t_mat(1,2) * imag_v(2,i)   &
                                         +               &
                       cosine_t_mat(1,3) * imag_v(3,i)   &
                                         +               &
                       cosine_t_mat(1,4) * imag_v(4,i)   &
                                         +               &
                       sine_t_mat(1,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(1,2)   * real_v(2,i)   &
                                         +               & 
                       sine_t_mat(1,3)   * real_v(3,i)   &
                                         +               &
                       sine_t_mat(1,4)   * real_v(4,i)  
     imag_v_scr(2,i) = cosine_t_mat(2,1) * imag_v(1,i)   &
                                         +               &
                       cosine_t_mat(2,2) * imag_v(2,i)   &
                                         +               &
                       cosine_t_mat(2,3) * imag_v(3,i)   &
                                         +               &
                       cosine_t_mat(2,4) * imag_v(4,i)   &
                                         +               &
                       sine_t_mat(2,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(2,2)   * real_v(2,i)   &
                                         +               &
                       sine_t_mat(2,3)   * real_v(3,i)   &
                                         +               &
                       sine_t_mat(2,4)   * real_v(4,i)  
     imag_v_scr(3,i) = cosine_t_mat(3,1) * imag_v(1,i)   &
                                         +               &
                       cosine_t_mat(3,2) * imag_v(2,i)   &
                                         +               &
                       cosine_t_mat(3,3) * imag_v(3,i)   &
                                         +               & 
                       cosine_t_mat(3,4) * imag_v(4,i)   &
                                         +               &
                       sine_t_mat(3,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(3,2)   * real_v(2,i)   &
                                         +               &
                       sine_t_mat(3,3)   * real_v(3,i)   &
                                         +               &
                       sine_t_mat(3,4)   * real_v(4,i)  
     imag_v_scr(4,i) = cosine_t_mat(4,1) * imag_v(1,i)   &
                                         +               &
                       cosine_t_mat(4,2) * imag_v(2,i)   &
                                         +               &
                       cosine_t_mat(4,3) * imag_v(3,i)   &
                                         +               &
                       cosine_t_mat(4,4) * imag_v(4,i)   &
                                         +               &
                       sine_t_mat(4,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(4,2)   * real_v(2,i)   &
                                         +               &
                       sine_t_mat(4,3)   * real_v(3,i)   &
                                         +               &
                       sine_t_mat(4,4)   * real_v(4,i)  
  END DO
  real_v(1:4,:) = real_v_scr(1:4,:)
  imag_v(1:4,:) = imag_v_scr(1:4,:)
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
  SUBROUTINE v_m_v_5(real_v,imag_v,real_v_scr,imag_v_scr, &
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
                                         -               &
                       sine_t_mat(1,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(1,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(1,3)   * imag_v(3,i)   &
                                         -               &
                       sine_t_mat(1,4)   * imag_v(4,i)   &
                                         -               &
                       sine_t_mat(1,5)   * imag_v(5,i)   
     real_v_scr(2,i) = cosine_t_mat(2,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(2,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(2,3) * real_v(3,i)   &
                                         +               &
                       cosine_t_mat(2,4) * real_v(4,i)   &
                                         +               &
                       cosine_t_mat(2,5) * real_v(5,i)   &
                                         -               &
                       sine_t_mat(2,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(2,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(2,3)   * imag_v(3,i)   &
                                         -               &
                       sine_t_mat(2,4)   * imag_v(4,i)   &
                                         -               &
                       sine_t_mat(2,5)   * imag_v(5,i)   
     real_v_scr(3,i) = cosine_t_mat(3,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(3,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(3,3) * real_v(3,i)   &
                                         +               &
                       cosine_t_mat(3,4) * real_v(4,i)   &
                                         +               &
                       cosine_t_mat(3,5) * real_v(5,i)   &
                                         -               &
                       sine_t_mat(3,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(3,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(3,3)   * imag_v(3,i)   &
                                         -               &
                       sine_t_mat(3,4)   * imag_v(4,i)   &
                                         -               &
                       sine_t_mat(3,5)   * imag_v(5,i)   
     real_v_scr(4,i) = cosine_t_mat(4,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(4,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(4,3) * real_v(3,i)   &
                                         +               &
                       cosine_t_mat(4,4) * real_v(4,i)   &
                                         +               &
                       cosine_t_mat(4,5) * real_v(5,i)   &  
                                         -               &
                       sine_t_mat(4,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(4,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(4,3)   * imag_v(3,i)   &
                                         -               &
                       sine_t_mat(4,4)   * imag_v(4,i)   &
                                         -               &
                       sine_t_mat(4,5)   * imag_v(5,i)   
     real_v_scr(5,i) = cosine_t_mat(5,1) * real_v(1,i)   &
                                         +               &
                       cosine_t_mat(5,2) * real_v(2,i)   &
                                         +               &
                       cosine_t_mat(5,3) * real_v(3,i)   &
                                         +               &
                       cosine_t_mat(5,4) * real_v(4,i)   &
                                         +               &
                       cosine_t_mat(5,5) * real_v(5,i)   &  
                                         -               &
                       sine_t_mat(5,1)   * imag_v(1,i)   &
                                         -               &
                       sine_t_mat(5,2)   * imag_v(2,i)   &
                                         -               &
                       sine_t_mat(5,3)   * imag_v(3,i)   &
                                         -               &
                       sine_t_mat(5,4)   * imag_v(4,i)   &
                                         -               &
                       sine_t_mat(5,5)   * imag_v(5,i)   
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
                       sine_t_mat(1,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(1,2)   * real_v(2,i)   &
                                         +               & 
                       sine_t_mat(1,3)   * real_v(3,i)   &
                                         +               &
                       sine_t_mat(1,4)   * real_v(4,i)   &
                                         +               &
                       sine_t_mat(1,5)   * real_v(5,i)   
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
                       sine_t_mat(2,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(2,2)   * real_v(2,i)   &
                                         +               &
                       sine_t_mat(2,3)   * real_v(3,i)   &
                                         +               &
                       sine_t_mat(2,4)   * real_v(4,i)   &
                                         +               &
                       sine_t_mat(2,5)   * real_v(5,i)  
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
                       sine_t_mat(3,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(3,2)   * real_v(2,i)   &
                                         +               &
                       sine_t_mat(3,3)   * real_v(3,i)   &
                                         +               &
                       sine_t_mat(3,4)   * real_v(4,i)   &
                                         +               &
                       sine_t_mat(3,5)   * real_v(5,i)  
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
                       sine_t_mat(4,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(4,2)   * real_v(2,i)   &
                                         +               &
                       sine_t_mat(4,3)   * real_v(3,i)   &
                                         +               &
                       sine_t_mat(4,4)   * real_v(4,i)   &
                                         +               &
                       sine_t_mat(4,5)   * real_v(5,i)  
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
                       sine_t_mat(5,1)   * real_v(1,i)   &
                                         +               &
                       sine_t_mat(5,2)   * real_v(2,i)   &
                                         +               &
                       sine_t_mat(5,3)   * real_v(3,i)   &
                                         +               &
                       sine_t_mat(5,4)   * real_v(4,i)   &
                                         +               &
                       sine_t_mat(5,5)   * real_v(5,i)  
  END DO
  real_v(1:5,:) = real_v_scr(1:5,:)
  imag_v(1:5,:) = imag_v_scr(1:5,:)
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
  SUBROUTINE v_m_v_6(real_v,imag_v,real_v_scr,imag_v_scr, &
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
                       sine_t_mat(1,6)   * imag_v(6,i)   
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
                       sine_t_mat(2,6)   * imag_v(6,i)   
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
                       sine_t_mat(3,6)   * imag_v(6,i)   
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
                       sine_t_mat(4,6)   * imag_v(6,i)   
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
                       sine_t_mat(5,6)   * imag_v(6,i)   
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
                       sine_t_mat(6,6)   * imag_v(6,i)   

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
                       sine_t_mat(1,6)   * real_v(6,i)   
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
                       sine_t_mat(2,6)   * real_v(6,i)   
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
                       sine_t_mat(3,6)   * real_v(6,i)   
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
                       sine_t_mat(4,6)   * real_v(6,i)   
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
                       sine_t_mat(5,6)   * real_v(6,i)   
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
                       sine_t_mat(6,6)   * real_v(6,i)   

  END DO
  real_v(1:6,:) = real_v_scr(1:6,:)
  imag_v(1:6,:) = imag_v_scr(1:6,:)
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
  SUBROUTINE v_m_v_8(real_v,imag_v,real_v_scr,imag_v_scr, &
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
                       sine_t_mat(1,8)   * imag_v(8,i)   
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
                       sine_t_mat(2,8)   * imag_v(8,i)   
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
                       sine_t_mat(3,8)   * imag_v(8,i)   
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
                       sine_t_mat(4,8)   * imag_v(8,i)   
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
                       sine_t_mat(5,8)   * imag_v(8,i)   
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
                       sine_t_mat(6,8)   * imag_v(8,i)   
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
                       sine_t_mat(7,8)   * imag_v(8,i)   
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
                       sine_t_mat(8,8)   * imag_v(8,i)   
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
                       sine_t_mat(1,8)   * real_v(8,i)   
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
                       sine_t_mat(2,8)   * real_v(8,i)   
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
                       sine_t_mat(3,8)   * real_v(8,i)   
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
                       sine_t_mat(4,8)   * real_v(8,i)   
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
                       sine_t_mat(5,8)   * real_v(8,i)   
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
                       sine_t_mat(6,8)   * real_v(8,i)   
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
                       sine_t_mat(7,8)   * real_v(8,i)   
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
                       sine_t_mat(8,8)   * real_v(8,i)   
  END DO
  real_v(1:8,:) = real_v_scr(1:8,:)
  imag_v(1:8,:) = imag_v_scr(1:8,:)
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
!
!
END MODULE v_m_v
