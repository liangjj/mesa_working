!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    MODULE finite_element_matrix_multiply
!deck finite_element_matrix_multiply
!***begin prologue     finite_element_matrix_multiply     
!***date written       040607   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            matrix multiply routine for DVR matrices.
!
!                     H = H(1,1) + H(2,2) + H(3,3)
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
!***                   
!***                   
!
!***end prologue       finite_element_matrix_multiply
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck dvr_m_v
!***begin prologue     dvr_m_v     
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
!                 V = M * V
!
!------------------------------------------------------------------------------------
!***                   Where M is a DVR matrix.
!
!***                   The parameter nj is a dummy and can take on values
!***                   consistent with any dimensional problem.  In a 1D
!***                   case nj=1, in 2D nj=nx and in 3D nj=ny*nx.
!***                   It is assumed that double counting of the diagonals
!***                   has been accounted for by halving the elements at
!***                   the interior, connecting, DVR element points.
!***references
!***routines called    
!
!***end prologue       dvr_m_v                     
!
  SUBROUTINE dvr_m_v(v,                    &
                     v_scr,                &
                     ni,nj,index)
  USE io
  USE dvr_global,                     ONLY    : num_reg, nfun_reg
  USE dvrprop_global_it,              ONLY    : mat_reg_d
  USE v_m_v
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, index
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  INTEGER                                  :: i
  INTEGER                                  :: locate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
!  The sequence is called three times.  First for the odd matrices, then
!  for the even matrices and then again for the odd matrices.
! 
  locate = 1
  DO i = 1, num_reg(index)
!
    IF ( nfun_reg(i,index) > 7 ) then
!
!                    General Code
!
         CALL v_m_v_gen(v(locate,1),                    &
                        v_scr(locate,1),                &
                        mat_reg_d(i,index)%ke_mat_d,    &
                        ni,nj,nfun_reg(i,index))
    ELSE IF (nfun_reg(i,index) == 2) then
! 
!                       Special case for Two by Two
!
        CALL v_m_v_2(v(locate,1),                       &
                     v_scr(locate,1),                   &
                     mat_reg_d(i,index)%ke_mat_d,       &
                     ni,nj)
    ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Special case for Three by Three
!
        CALL v_m_v_3(v(locate,1),                       &
                     v_scr(locate,1),                   &
                     mat_reg_d(i,index)%ke_mat_d,       &
                     ni,nj)
    ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Special case for Four by Four
!
        CALL v_m_v_4(v(locate,1),                       &
                     v_scr(locate,1),                   &
                     mat_reg_d(i,index)%ke_mat_d,       &
                     ni,nj)
    ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Special case for Five by Five
!
        CALL v_m_v_5(v(locate,1),                       &
                     v_scr(locate,1),                   &
                     mat_reg_d(i,index)%ke_mat_d,       &
                     ni,nj)
    ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Special case for Six by Six
!
        CALL v_m_v_6(v(locate,1),                      &
                     v_scr(locate,1),                  &
                     mat_reg_d(i,index)%ke_mat_d,      &
                     ni,nj)
    ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Special case for Seven by Seven
!
        CALL v_m_v_7(v(locate,1),                      &
                     v_scr(locate,1),                  &
                     mat_reg_d(i,index)%ke_mat_d,      &
                     ni,nj)
    ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Special case for Eight by Eight
!
        CALL v_m_v_8(v(locate,1),                      &
                     v_scr(locate,1),                  &
                     mat_reg_d(i,index)%ke_mat_d,      &
                        ni,nj)
    ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Special case for Nine by Nine
!
        CALL v_m_v_9(v(locate,1),                      &
                     v_scr(locate,1),                  &
                     mat_reg_d(i,index)%ke_mat_d,      &
                     ni,nj)
    ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Special case for Ten by Ten
!
        CALL v_m_v_10(v(locate,1),                     &
                      v_scr(locate,1),                 &
                      mat_reg_d(i,index)%ke_mat_d,     &
                      ni,nj)
    END IF
    locate = locate + nfun_reg(i,index)  - 1
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE dvr_m_v
!
!

!*deck dvr_v_m_2_d
!***begin prologue     dvr_v_m_2_d     
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
!                 V = V * M 
!
!------------------------------------------------------------------------------------
!***                   where M is the regional
!***                   propagator for a FEDVR or FD Hamiltonian. 
!***                   parts.  
!***                   The routine is not needed in 1D problems. 
!***                   In 2D, the ni index runs over the number
!                      of points in the y coordinate while in 3D,
!***                   it runs over the product nz*ny of the number 
!***                   of points in the z and y coordinates.
!
!***references
!***routines called    
!***end prologue       dvr_v_m_2_d
!
  SUBROUTINE dvr_v_m_2_d(v,                    &
                         v_scr,                &
                         ni,nj,index)
  USE io
  USE dvr_global,                ONLY    : num_reg, nfun_reg
  USE dvrprop_global_it,         ONLY    : mat_reg_d
  USE v_v_m
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, index
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  INTEGER                                  :: i
  INTEGER                                  :: locate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  locate = 1
!
  DO i=1, num_reg(index)
     IF ( nfun_reg(i,index) > 7 ) then
!
!                        General Case
!
          CALL v_v_m_gen(v(1,locate),                     &
                         v_scr(1,locate),                 &
                         mat_reg_d(i,index)%ke_mat_d,     &
                         ni,nj,nfun_reg(i,index)) 
!
     ELSE IF (nfun_reg(i,index) == 2) then
!
!                       Two by Two
! 
          CALL v_v_m_2(v(1,locate),                       &
                       v_scr(1,locate),                   &
                       mat_reg_d(i,index)%ke_mat_d,       &
                       ni,nj)
     ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Three by Three
!
          CALL v_v_m_3(v(1,locate),                       &
                       v_scr(1,locate),                   &
                       mat_reg_d(i,index)%ke_mat_d,       &
                       ni,nj)
     ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Four by Four
!
          CALL v_v_m_4(v(1,locate),                       &
                       v_scr(1,locate),                   &
                       mat_reg_d(i,index)%ke_mat_d,       &
                       ni,nj)
     ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Five by Five
!
          CALL v_v_m_5(v(1,locate),                       &
                       v_scr(1,locate),                   &
                       mat_reg_d(i,index)%ke_mat_d,       &
                       ni,nj)
     ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Six by Six
!
          CALL v_v_m_6(v(1,locate),                       &
                       v_scr(1,locate),                   &
                       mat_reg_d(i,index)%ke_mat_d,       &
                       ni,nj)
     ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Seven by Seven
!
          CALL v_v_m_7(v(1,locate),                       &
                       v_scr(1,locate),                   &
                       mat_reg_d(i,index)%ke_mat_d,       &
                       ni,nj)
     ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Eight by Eight
!
          CALL v_v_m_8(v(1,locate),                       &
                       v_scr(1,locate),                   &
                       mat_reg_d(i,index)%ke_mat_d,       &
                       ni,nj)
     ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Nine by Nine
!
          CALL v_v_m_9(v(1,locate),                       &
                       v_scr(1,locate),                   &
                       mat_reg_d(i,index)%ke_mat_d,       &
                       ni,nj)
     ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Ten by Ten
!
          CALL v_v_m_10(v(1,locate),                      &
                        v_scr(1,locate),                  &
                        mat_reg_d(i,index)%ke_mat_d,      &
                        ni,nj)
     END IF
     locate = locate + nfun_reg(i,index)  - 1
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE dvr_v_m_2_d
!
!
!*deck dvr_v_m_3_d
!***begin prologue     dvr_v_m_3_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine simply calls
!***                   dvr_v_m_2_d in a loop over a dummy
!***                   index, k
!
!***references
!***routines called
!***end prologue       dvr_v_m_3_d
!
  SUBROUTINE dvr_v_m_3_d(v,                   &
                         v_scr,               &
                         ni,nj,nk,index)
  IMPLICIT NONE
  INTEGER                             :: ni, nj, nk, index
  REAL*8, DIMENSION(ni,nj,nk)         :: v
  REAL*8, DIMENSION(ni,nj,nk)         :: v_scr
  INTEGER                             :: k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  DO k=1,nk
     CALL dvr_v_m_2_d(v(1,1,k),                  &
                      v_scr(1,1,k),              &
                      ni,nj,index) 
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE dvr_v_m_3_d
END MODULE finite_element_matrix_multiply

