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
