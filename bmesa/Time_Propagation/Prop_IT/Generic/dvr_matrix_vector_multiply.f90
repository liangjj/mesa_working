!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     MODULE dvr_matrix_vector_multiply
                        USE io
                        USE dvr_global,         ONLY   : num_reg, nfun_reg
                        USE v_out_eq_dvr_mat_v_in
                        USE v_out_eq_v_in_dvr_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck dvr_mat_mul_d
!***begin prologue     dvr_mat_mul_d    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine computes the matrix vector multiply 
!***                   needed for a Hamiltonian which is of the FEDVR
!**                    form.  There is no need for a packed matrix here
!**                    as the structure of the M matrix is recognized.
!**                    The double counting of the interface elements must
!**                    taken care of before calling this routine.  Other
!**                    routines are used to do this by halving the diagonal
!**                    elements to avoid overcounting. 
!------------------------------------------------------------------------------------
!
!                             V = M * V
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
!**                    Special routines have been written for matrices
!**                    up to 10*10 for efficiency.
!***references
!***routines called    
!
!***end prologue       dvr_mat_mul_d                     
!
  SUBROUTINE dvr_mat_mul_d(v,v_scr,ni,nj,index)
  USE dvrprop_global_it,  ONLY   : mat_reg_d
  IMPLICIT NONE
  INTEGER                                :: ni, nj, index
  REAL*8, DIMENSION(ni,nj)               :: v
  REAL*8, DIMENSION(ni,nj)               :: v_scr
  INTEGER                                :: i
  INTEGER                                :: locate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  locate = 1
  DO i = 1, num_reg(index)
!
    IF ( nfun_reg(i,index) > 10 ) then
!
!                    General Code
!
         CALL v_mat_v_gen(v(locate:,:),                  &
                          v_scr(locate:,:),              &
                          mat_reg_d(i,index)%ke_mat_d,   &
                          ni,nj,nfun_reg(i,index))
    ELSE IF (nfun_reg(i,index) == 2) then
! 
!                       Special case for Two by Two
!
        CALL v_mat_v_2(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg_d(i,index)%ke_mat_d,      &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Special case for Three by Three
!
        CALL v_mat_v_3(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg_d(i,index)%ke_mat_d,      &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Special case for Four by Four
!
        CALL v_mat_v_4(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg_d(i,index)%ke_mat_d,      &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Special case for Five by Five
!
        CALL v_mat_v_5(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg_d(i,index)%ke_mat_d,      &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Special case for Six by Six
!
        CALL v_mat_v_6(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg_d(i,index)%ke_mat_d,      &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Special case for Seven by Seven
!
        CALL v_mat_v_7(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg_d(i,index)%ke_mat_d,      &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Special case for Eight by Eight
!
        CALL v_mat_v_8(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg_d(i,index)%ke_mat_d,      &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Special case for Nine by Nine
!
        CALL v_mat_v_9(v(locate:,:),                     &
                       v_scr(locate:,:),                 &
                       mat_reg_d(i,index)%ke_mat_d,      &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Special case for Ten by Ten
!
        CALL v_mat_v_10(v(locate:,:),                    &
                        v_scr(locate:,:),                &
                        mat_reg_d(i,index)%ke_mat_d,     &
                        ni,nj)
    END IF
    locate = locate + nfun_reg(i,index)  - 1
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE dvr_mat_mul_d
!
!deck dvr_mat_mul_z
!***begin prologue     dvr_mat_mul_z    
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
!                             V = M * V
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
!***end prologue       dvr_mat_mul_z                     
!
  SUBROUTINE dvr_mat_mul_z(v,v_scr,ni,nj,index)
  USE dvrprop_global_rt,  ONLY   : mat_reg_d
  IMPLICIT NONE
  INTEGER                                :: ni, nj, index
  COMPLEX*16, DIMENSION(ni,nj)           :: v
  COMPLEX*16, DIMENSION(ni,nj)           :: v_scr
  INTEGER                                :: i
  INTEGER                                :: locate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  locate = 1
  DO i = 1, num_reg(index)
!
    IF ( nfun_reg(i,index) > 10 ) then
!
!                    General Code
!
         CALL v_mat_v_gen(v(locate:,:),                         &
                          v_scr(locate:,:),                     &
                          mat_reg_d(i,index)%ke_mat_d,          &
                          ni,nj,nfun_reg(i,index))
    ELSE IF (nfun_reg(i,index) == 2) then
! 
!                       Special case for Two by Two
!
        CALL v_mat_v_2(v(locate:,:),                            &
                       v_scr(locate:,:),                        &
                       mat_reg_d(i,index)%ke_mat_d,             &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Special case for Three by Three
!
        CALL v_mat_v_3(v(locate:,:),                            &
                       v_scr(locate:,:),                        &
                       mat_reg_d(i,index)%ke_mat_d,             &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Special case for Four by Four
!
        CALL v_mat_v_4(v(locate:,:),                            &
                       v_scr(locate:,:),                        &
                       mat_reg_d(i,index)%ke_mat_d,             &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Special case for Five by Five
!
        CALL v_mat_v_5(v(locate:,:),                            &
                       v_scr(locate:,:),                        &
                       mat_reg_d(i,index)%ke_mat_d,             &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Special case for Six by Six
!
        CALL v_mat_v_6(v(locate:,:),                            &
                       v_scr(locate:,:),                        &
                       mat_reg_d(i,index)%ke_mat_d,             &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Special case for Seven by Seven
!
        CALL v_mat_v_7(v(locate:,:),                            &
                       v_scr(locate:,:),                        &
                       mat_reg_d(i,index)%ke_mat_d,             &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Special case for Eight by Eight
!
        CALL v_mat_v_8(v(locate:,:),                            &
                       v_scr(locate:,:),                        &
                       mat_reg_d(i,index)%ke_mat_d,             &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Special case for Nine by Nine
!
        CALL v_mat_v_9(v(locate:,:),                            &
                       v_scr(locate:,:),                        &
                       mat_reg_d(i,index)%ke_mat_d,             &
                       ni,nj)
    ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Special case for Ten by Ten
!
        CALL v_mat_v_10(v(locate:,:),                           &
                        v_scr(locate:,:),                       &
                        mat_reg_d(i,index)%ke_mat_d,            &
                        ni,nj)
    END IF
    locate = locate + nfun_reg(i,index)  - 1
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE dvr_mat_mul_z
!
!*deck dvr_mat_mul_2_d
!***begin prologue     dvr_mat_mul_2_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine is the driver routine which
!***                   calculates the full two dimensional matrix
!**                    vector multiply using the fact that the matrices
!**                    are separable in each dimension.
!------------------------------------------------------------------------------------
!
!                            V = V * M 
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
!***end prologue       dvr_mat_mul_2_d
!
  SUBROUTINE dvr_mat_mul_2_d(v,v_scr,ni,nj,nv,index)
  USE dvrprop_global_it,  ONLY   : mat_reg_d
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, nv, index
  REAL*8, DIMENSION(ni,nj,nv)              :: v
  REAL*8, DIMENSION(ni,nj,nv)              :: v_scr
  INTEGER                                  :: i, j
  INTEGER                                  :: locate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  locate = 1
!
  DO i=1, num_reg(index)
     IF ( nfun_reg(i,index) > 10 ) then
!
!                        General Case
!
          DO j=1,nv
             CALL v_v_mat_gen(v(:,locate:,j),                    &
                               v_scr(:,locate:,j),               &
                               mat_reg_d(i,index)%ke_mat_d,      &
                               ni,nj,nfun_reg(i,index)) 
          END DO
!
     ELSE IF (nfun_reg(i,index) == 2) then
!
!                       Two by Two
! 
          DO j=1,nv
             CALL v_v_mat_2(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg_d(i,index)%ke_mat_d,        &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Three by Three
!
          DO j=1,nv
             CALL v_v_mat_3(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg_d(i,index)%ke_mat_d,        &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Four by Four
!
          DO j=1,nv
             CALL v_v_mat_4(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg_d(i,index)%ke_mat_d,        &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Five by Five
!
          DO j=1,nv
             CALL v_v_mat_5(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg_d(i,index)%ke_mat_d,        &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Six by Six
!
          DO j=1,nv
             CALL v_v_mat_6(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg_d(i,index)%ke_mat_d,        &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Seven by Seven
!
          DO j=1,nv
             CALL v_v_mat_7(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg_d(i,index)%ke_mat_d,        &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Eight by Eight
!
          DO j=1,nv
             CALL v_v_mat_8(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg_d(i,index)%ke_mat_d,        &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Nine by Nine
!
          DO j=1,nv
             CALL v_v_mat_9(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg_d(i,index)%ke_mat_d,        &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Ten by Ten
!
          DO j=1,nv
             CALL v_v_mat_10(v(:,locate:,j),                    &
                             v_scr(:,locate:,j),                &
                             mat_reg_d(i,index)%ke_mat_d,       &
                             ni,nj)
          END DO
     END IF
     locate = locate + nfun_reg(i,index)  - 1
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE dvr_mat_mul_2_d
!
!*deck dvr_mat_mul_2_z
!***begin prologue     dvr_mat_mul_2_z     
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
!***end prologue       dvr_mat_mul_2_z
!
  SUBROUTINE dvr_mat_mul_2_z(v,v_scr,ni,nj,nv,index)
  USE dvrprop_global_rt,  ONLY   : mat_reg_d
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, nv, index
  COMPLEX*16, DIMENSION(ni,nj,nv)          :: v
  COMPLEX*16, DIMENSION(ni,nj,nv)          :: v_scr
  INTEGER                                  :: i, j
  INTEGER                                  :: locate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  locate = 1
!
  DO i=1, num_reg(index)
     IF ( nfun_reg(i,index) > 10 ) then
!
!                        General Case
!
          DO j=1,nv
             CALL v_v_mat_gen(v(:,locate:,j),                   &
                              v_scr(:,locate:,j),               &
                              mat_reg_d(i,index)%ke_mat_d,      &
                              ni,nj,nfun_reg(i,index)) 
          END DO
!
     ELSE IF (nfun_reg(i,index) == 2) then
!
!                       Two by Two
! 
          DO j=1,nv
             CALL v_v_mat_2(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg_d(i,index)%ke_mat_d,        &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Three by Three
!
          DO j=1,nv
             CALL v_v_mat_3(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg_d(i,index)%ke_mat_d,        &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Four by Four
!
          DO j=1,nv
             CALL v_v_mat_4(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg_d(i,index)%ke_mat_d,        &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Five by Five
!
          DO j=1,nv
             CALL v_v_mat_5(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg_d(i,index)%ke_mat_d,        &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Six by Six
!
          DO j=1,nv
             CALL v_v_mat_6(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg_d(i,index)%ke_mat_d,        &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Seven by Seven
!
          DO j=1,nv
             CALL v_v_mat_7(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg_d(i,index)%ke_mat_d,        &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Eight by Eight
!
          DO j=1,nv
             CALL v_v_mat_8(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg_d(i,index)%ke_mat_d,        &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Nine by Nine
!
          DO j=1,nv
             CALL v_v_mat_9(v(:,locate:,j),                     &
                            v_scr(:,locate:,j),                 &
                            mat_reg_d(i,index)%ke_mat_d,        &
                            ni,nj)
          END DO
     ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Ten by Ten
!
          DO j=1,nv
             CALL v_v_mat_10(v(:,locate:,j),                    &
                             v_scr(:,locate:,j),                &
                             mat_reg_d(i,index)%ke_mat_d,       &
                             ni,nj)
          END DO
     END IF
     locate = locate + nfun_reg(i,index)  - 1
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE dvr_mat_mul_2_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               END   MODULE dvr_matrix_vector_multiply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
