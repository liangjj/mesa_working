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
