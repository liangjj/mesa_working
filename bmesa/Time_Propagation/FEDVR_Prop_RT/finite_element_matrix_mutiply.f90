!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   MODULE finite_element_matrix_multiply
!deck finite_element_matrix_multiply
!***begin prologue     finite_element_matrix_multiply     
!***date written       040607   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            compute the effect of the split operator
!                      off-diagonal, exponential propagator on a vector.
!
!***description        The operation consists of applying a DVR matrix
!***                   to a vector.
!
!                      
!                      H = H(1,1) + H(2,2) + H(3,3)
!
!     Note that vectors are stored as V(nphy(3),nphy(2),nphy(1)).  Thus matrix
!     vector multiply in two and three dimensions must be done carefully for
!     nphy(2) and nphy(3).  We have( repeated indices summed over), in 3D,
!    
! V(3,2,1) = [ H(1,1) * V(3,2,1) + H(2,2) * V(3,2,1) + H(3,3) * V(3,2,1) ]  
!                                   =
!            [ V(3,2,1) * H(1,1) + V(3,2,1) * H(2,2) + H(3,3) * V(3,2,1) ]
!
! Thus the first mutiply may be done as a matrix V(3*2,1) on H(1,1), 
! the second as V(3,2,1) * H(2,2) with an outer loop over index 1 and 
! the third as a simple matrix vector multiply, H(3,3) * V(3,2*1).
!
!and in 2D,
!
! V(2,1)   = [ H(1,1) * V(2,1) + H(2,2) * V(2,1)  ]  
!                                   =
!            [ V(2,1) * H(1,1) + H(2,2) * V(2,1) ]
!
!
!
!***references
!***routines called    
!
!***end prologue       finite_element_matrix_multiply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        INTERFACE dvr_mat_mul
             MODULE PROCEDURE dvr_mat_mul_1d_z, dvr_mat_mul_2d,z    &
                              dvr_mat_mul_3d_z
                    END INTERFACE dvr_mat_mul
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!deck dvr_mat_mul_m_v
!***begin prologue     dvr_mat_mul_m_v     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine is the driver routine which 
!***                   calculates,
!------------------------------------------------------------------------------------
!
!                 Real(V) = M * Real(V) 
!                 Imag(V) = M * Imag(V) 
!
!------------------------------------------------------------------------------------
!***                   The parameter nj is a dummy and can take on values
!***                   consistent with any dimensional problem.  In a 1D
!***                   case nj=1, in 2D nj=nx and in 3D nj=ny*nx.

!***references
!***routines called    
!
!***end prologue       dvr_mat_mul_m_v                     
!
  SUBROUTINE dvr_mat_mul_m_v(real_v,             &
                             imag_v,             &
                             real_v_scr,         &
                             imag_v_scr,         &
                             ni,nj,index)
  USE io
  USE dvr_global,                  ONLY    : num_reg, nfun_reg
  USE dvrprop_global_rt,           ONLY    : mat_reg_d, prop_point
  USE v_m_v
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, index
  REAL*8, DIMENSION(ni,nj)                 :: real_v, imag_v
  REAL*8, DIMENSION(ni,nj)                 :: real_v_scr, imag_v_scr
  INTEGER                                  :: i, ntrips, trips
  INTEGER, DIMENSION(3)                    :: locate, begin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  locate = 1
  DO i = 1, num_reg(index)
!
     IF ( nfun_reg(i,index) > 7 ) then
!
!                       General Code
!
         CALL v_m_v_gen(real_v(locate(trips),1),                              &
                          imag_v(locate(trips),1),                              &
                          real_v_scr(locate(trips),1),                          &
                          imag_v_scr(locate(trips),1),                          &
                          mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),      &
                          mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),        &
                          ni,nj,nfun_reg(i,index))
        ELSE IF (nfun_reg(i,index) == 2) then
! 
!                       Special case for Two by Two
!
           CALL v_m_v_2(real_v(locate(trips),1),                                &
                        imag_v(locate(trips),1),                                &
                        real_v_scr(locate(trips),1),                            &
                        imag_v_scr(locate(trips),1),                            &
                        mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),        &
                        mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),          &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Special case for Three by Three
!
           CALL v_m_v_3(real_v(locate(trips),1),                                &
                        imag_v(locate(trips),1),                                &
                        real_v_scr(locate(trips),1),                            &
                        imag_v_scr(locate(trips),1),                            &
                        mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),        &
                        mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),          &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Special case for Four by Four
!
           CALL v_m_v_4(real_v(locate(trips),1),                                &
                        imag_v(locate(trips),1),                                &
                        real_v_scr(locate(trips),1),                            &
                        imag_v_scr(locate(trips),1),                            &
                        mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),        &
                        mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),          &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Special case for Five by Five
!
           CALL v_m_v_5(real_v(locate(trips),1),                                &
                        imag_v(locate(trips),1),                                &
                        real_v_scr(locate(trips),1),                            &
                        imag_v_scr(locate(trips),1),                            &
                        mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),        &
                        mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),          &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Special case for Six by Six
!
           CALL v_m_v_6(real_v(locate(trips),1),                                &
                        imag_v(locate(trips),1),                                &
                        real_v_scr(locate(trips),1),                            &
                        imag_v_scr(locate(trips),1),                            &
                        mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),        &
                        mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),          &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Special case for Seven by Seven
!
           CALL v_m_v_7(real_v(locate(trips),1),                                &
                        imag_v(locate(trips),1),                                &
                        real_v_scr(locate(trips),1),                            &
                        imag_v_scr(locate(trips),1),                            &
                        mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),        &
                        mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),          &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Special case for Eight by Eight
!
           CALL v_m_v_8(real_v(locate(trips),1),                                &
                        imag_v(locate(trips),1),                                &
                        real_v_scr(locate(trips),1),                            &
                        imag_v_scr(locate(trips),1),                            &
                        mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),        &
                        mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),          &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Special case for Nine by Nine
!
           CALL v_m_v_9(real_v(locate(trips),1),                                &
                        imag_v(locate(trips),1),                                &
                        real_v_scr(locate(trips),1),                            &
                        imag_v_scr(locate(trips),1),                            &
                        mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),        &
                        mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),          &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Special case for Ten by Ten
!
           CALL v_m_v_10(real_v(locate(trips),1),                               &
                         imag_v(locate(trips),1),                               &
                         real_v_scr(locate(trips),1),                           &
                         imag_v_scr(locate(trips),1),                           &
                         mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),       &
                         mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),         &
                         ni,nj)
        END IF
        IF (i /= num_reg(index) ) then
            locate(trips) = locate(trips) + nfun_reg(i,index)   &
                                          +                     &
                            nfun_reg(i+1,index)                 &
                                          -                     &
                                          2
        END IF
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE dvr_mat_mul_m_v
!
!

!*deck dvr_mat_mul_v_m_2_d
!***begin prologue     dvr_mat_mul_v_m_2_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine is the driver routine which
!***                   calculates,
!------------------------------------------------------------------------------------
!
!                 Real(V) = Real(V) * cos_t_mat - Imag(V) * sin_t_mat
!                 Imag(V) = Imag(V) * cos_t_mat + Real(V) * sin_t_mat
!
!------------------------------------------------------------------------------------
!***                   where cos_t_mat and sine_t_mat_2 are the regional
!***                   propagators for a FEDVR or FD Hamiltonian. 
!***                   parts.  
!***                   The routine is not needed in 1D problems. 
!***                   In 2D, the ni index runs over the number
!                      of points in the y coordinate while in 3D,
!***                   it runs over the product nz*ny of the number 
!***                   of points in the z and y coordinates.
!
!***references
!***routines called    v_v_m_gen, v_v_m_2, v_v_m_3
!***end prologue       dvr_mat_mul_v_m_2_d
!
  SUBROUTINE dvr_mat_mul_v_m_2_d(real_v,       &
                                  imag_v,       &
                                  real_v_scr,   &
                                  imag_v_scr,   &
                                  ni,nj,index)
  USE io
  USE dvr_global,             ONLY    : num_reg, nfun_reg
  USE dvrprop_global_rt,      ONLY    : mat_reg_d, prop_point
  USE v_v_m
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, index, ntrips
  REAL*8, DIMENSION(ni,nj)                 :: real_v, imag_v
  REAL*8, DIMENSION(ni,nj)                 :: real_v_scr, imag_v_scr
  INTEGER                                  :: i, trips
  INTEGER, DIMENSION(3)                    :: locate, begin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  locate(1) = 1
  locate(2) = nfun_reg(1,index)
  locate(3) = 1
  begin(1) = 1
  begin(2) = 2
  begin(3) = 1
  ntrips = 3
  IF(num_reg(index) == 1 ) THEN
     ntrips=1
  END IF
!
  DO trips = 1, ntrips
     DO i=begin(trips), num_reg(index), 2
        IF ( nfun_reg(i,index) > 7 ) then
!
!                        General Case
!
           CALL v_v_m_gen(real_v(1,locate(trips)),                               &
                          imag_v(1,locate(trips)),                               &
                          real_v_scr(1,locate(trips)),                           &
                          imag_v_scr(1,locate(trips)),                           &
                          mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),       &
                          mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),         &
                          ni,nj,nfun_reg(i,index)) 
!
        ELSE IF (nfun_reg(i,index) == 2) then
!
!                       Two by Two
! 
           CALL v_v_m_2(real_v(1,locate(trips)),                                 &
                        imag_v(1,locate(trips)),                                 &
                        real_v_scr(1,locate(trips)),                             &
                        imag_v_scr(1,locate(trips)),                             &
                        mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),         &
                        mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),           &
                        ni,nj)
       ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Three by Three
!
           CALL v_v_m_3(real_v(1,locate(trips)),                                 &
                        imag_v(1,locate(trips)),                                 &
                        real_v_scr(1,locate(trips)),                             &
                        imag_v_scr(1,locate(trips)),                             &   
                        mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),         &
                        mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),           &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Four by Four
!
           CALL v_v_m_4(real_v(1,locate(trips)),                                 &
                        imag_v(1,locate(trips)),                                 &
                        real_v_scr(1,locate(trips)),                             &
                        imag_v_scr(1,locate(trips)),                             &   
                        mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),         &
                        mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),           &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Five by Five
!
           CALL v_v_m_5(real_v(1,locate(trips)),                                 &
                        imag_v(1,locate(trips)),                                 &
                        real_v_scr(1,locate(trips)),                             &
                        imag_v_scr(1,locate(trips)),                             &   
                        mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),         &
                        mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),           &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Six by Six
!
           CALL v_v_m_6(real_v(1,locate(trips)),                                 &
                        imag_v(1,locate(trips)),                                 &
                        real_v_scr(1,locate(trips)),                             &
                        imag_v_scr(1,locate(trips)),                             &   
                        mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),         &
                        mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),           &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Seven by Seven
!
           CALL v_v_m_7(real_v(1,locate(trips)),                                 &
                        imag_v(1,locate(trips)),                                 &
                        real_v_scr(1,locate(trips)),                             &
                        imag_v_scr(1,locate(trips)),                             &   
                        mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),         &
                        mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),           &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Eight by Eight
!
           CALL v_v_m_8(real_v(1,locate(trips)),                                 &
                        imag_v(1,locate(trips)),                                 &
                        real_v_scr(1,locate(trips)),                             &
                        imag_v_scr(1,locate(trips)),                             &   
                        mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),         &
                        mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),           &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Nine by Nine
!
           CALL v_v_m_9(real_v(1,locate(trips)),                                 &
                        imag_v(1,locate(trips)),                                 &
                        real_v_scr(1,locate(trips)),                             &
                        imag_v_scr(1,locate(trips)),                             &   
                        mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),         &
                        mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),           &
                        ni,nj)
        ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Ten by Ten
!
           CALL v_v_m_10(real_v(1,locate(trips)),                                &
                         imag_v(1,locate(trips)),                                &
                         real_v_scr(1,locate(trips)),                            &
                         imag_v_scr(1,locate(trips)),                            &   
                         mat_reg_d(i,index)%cosine_t_mat(:,:,prop_point),        &
                         mat_reg_d(i,index)%sine_t_mat(:,:,prop_point),          &
                         ni,nj)
        END IF
        IF (i /= num_reg(index) ) then
            locate(trips) = locate(trips) + nfun_reg(i,index)                    &
                                          +                                      &
                            nfun_reg(i+1,index)                                  &
                                          -                                      &
                                          2
        END IF
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE dvr_mat_mul_v_m_2_d
!
!
!*deck dvr_mat_mul_v_m_3_d
!***begin prologue     dvr_mat_mul_v_m_3_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine simply calls
!***                   dvr_mat_mul_v_m_2_d in a loop over a dummy
!***                   index, k
!
!***references
!***routines called
!***end prologue       dvr_mat_mul_v_m_3_d
!
  SUBROUTINE dvr_mat_mul_v_m_3_d(real_v,     &
                                 imag_v,      &
                                 real_v_scr,  &
                                 imag_v_scr,  &
                                 ni,nj,nk,index)
  IMPLICIT NONE
  INTEGER                             :: ni, nj, nk, index
  REAL*8, DIMENSION(ni,nj,nk)         :: real_v, imag_v
  REAL*8, DIMENSION(ni,nj,nk)         :: real_v_scr, imag_v_scr
  INTEGER                             :: k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  DO k=1,nk
     CALL dvr_mat_mul_v_m_2_d(real_v(1,1,k),          &
                               imag_v(1,1,k),          &
                               real_v_scr(1,1,k),      &
                               imag_v_scr(1,1,k),      &
                               ni,nj,index) 
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE dvr_mat_mul_v_m_3_d
END MODULE finite_element_matrix_multiply

