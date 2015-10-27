!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dvr_matrix_vector_multiply_module
!**begin prologue     dvr_matrix_vector_multiply_module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Contains the routines needed to multiply the
!***                  Hamiltonian on a vector or to mutiply the exponential
!***                  sector propagators on a vector. Explicit interfaces are 
!***                  used to allow a transparent use of generic subroutines 
!***                  which work for both real and complex vectors.  
!***                  This feature permits a single code to be used for both 
!***                  real and imaginary time propagation.
!***description       See the specific routined.
!**references
!**modules needed     See USE statements below
!**end prologue       dvr_matrix_vector_multiply_module
!***********************************************************************
!***********************************************************************
                      MODULE dvr_matrix_vector_multiply_module
                         USE dvrprop_global
                         USE dvr_shared
                         USE dvr_global
                         USE Iterative_Global
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE finite_element_m_v
              MODULE PROCEDURE finite_element_m_v_d,                  &
                               finite_element_m_v_z
                END INTERFACE finite_element_m_v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE dvr_mat_mul
              MODULE PROCEDURE dvr_mat_mul_d,                           &
                               dvr_mat_mul_z
                END INTERFACE dvr_mat_mul
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE v_mat_v_gen
              MODULE PROCEDURE v_mat_v_gen_d,                           &
                               v_mat_v_gen_z
                END INTERFACE  v_mat_v_gen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE v_v_mat_gen
              MODULE PROCEDURE v_v_mat_gen_d,                           &
                               v_v_mat_gen_z
                END INTERFACE  v_v_mat_gen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE v_so_mat_v_gen
              MODULE PROCEDURE v_so_mat_v_gen_d,                        &
                               v_so_mat_v_gen_z
                END INTERFACE  v_so_mat_v_gen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE v_v_so_mat_gen
              MODULE PROCEDURE v_v_so_mat_gen_d,                        &
                               v_v_so_mat_gen_z
                END INTERFACE  v_v_so_mat_gen
                    INTERFACE v_v_so_mat_2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck finite_element_matrix_multiply
!***begin prologue     finite_element_matrix_multiply     
!***date written       040607   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Driving routine to matrix multiply routine for DVR matrices.
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck finite_element_m_v_d
!***begin prologue     finite_element_m_v_d    
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
!***end prologue       finite_element_m_v_d                     
!
  SUBROUTINE finite_element_m_v_d(v_in,v_out)
  IMPLICIT NONE
  INTEGER                                 :: i
  REAL*8, DIMENSION(:)                    :: v_in
  REAL*8, DIMENSION(:)                    :: v_out
!
!
  v_out(:) = 0.d0
  IF(spdim == 1) THEN
     CALL dvr_mat_mul_d(v_in,v_out,nphy(1),1,1)
     CALL v_diag_mul_1d_d(v_in,v_out)
  ELSE IF(spdim == 2) THEN
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       This code is for the general problem involving nphy(2),
!       that is H(2,2) * V(2,1) 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL dvr_mat_mul_d(v_in,v_out,nphy(2),nphy(1),2)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       This code is for the general problem involving nphy(1),
!       V(2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        CALL dvr_mat_mul_2d_d(v_in,v_out,nphy(2),nphy(1),1,1)
        CALL v_diag_mul_2d_d(v_in,v_out)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ELSE IF(spdim == 3) THEN
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            This code is for the general problem involving nphy(3),
!            H(3,3) * V(3,2,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        CALL dvr_mat_mul_d(v_in,v_out,nphy(3),nphy(2)*nphy(1),3)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general case involving the matrix multiply 
!        involving nphy(2), that is, V(3,2,1) * H(2,2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL dvr_mat_mul_2d_d(v_in,v_out,nphy(3),nphy(2),nphy(1),2)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general problem involving the matrix multiply 
!        for nphy(1), that is, V(3,2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL dvr_mat_mul_2d_d(v_in,v_out,nphy(3)*nphy(2),nphy(1),1,1)
        CALL v_diag_mul_3d_d(v_in,v_out)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END IF
  END SUBROUTINE finite_element_m_v_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
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
!**                    as the structure of the H matrix is recognized.
!**                    To avoid double counting of the interface elements 
!**                    which are the diagonal elements, we remove them
!**                    and then add them in after the rest of the multiplicatio
!**                    is finished.  In one dimension you need to think of this
!**                    as a set of overlapping matrices where the interface is the
!**                    the diagonal.  So, a simple example on a (4*4)would be;
!**
!**                                  XXX0
!**                                  XXX0 
!**                                  XXXX
!**                                  00XX
!**                    After zeroing the diagonal, the matrix-vector multiply may be performed
!**                    as two completely separate dense matrix multiplies.  The diagonal
!**                    is then trivially restored at the end.
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
  USE dvrprop_global,  ONLY   : mat_reg
  IMPLICIT NONE
  INTEGER                                :: ni, nj, index
  REAL*8, DIMENSION(ni,nj)               :: v
  REAL*8, DIMENSION(ni,nj)               :: v_scr
  INTEGER                                :: i, j
  INTEGER                                :: first, last
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  first = 1
  DO i = 1, num_reg(index)
!
    IF ( nfun_reg(i,index) > 10 ) then
!
!                    General Code
!        
         last = first + nfun_reg(i,index) - 1
         CALL v_mat_v_gen(v(first:last,1:nj),                  &
                          v_scr(first:last,1:nj),              &
                          mat_reg(i,index)%ke_mat_d)
    ELSE IF (nfun_reg(i,index) == 2) then
! 
!                       Special case for Two by Two
!
         last = first + 1
         CALL v_mat_v_2(v(first:last,1:nj),                     &
                        v_scr(first:last,1:nj),                 &
                        mat_reg(i,index)%ke_mat_d)

    ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Special case for Three by Three
         last = first + 2
         CALL v_mat_v_3(v(first:last,1:nj),                     &
                       v_scr(first:last,1:nj),                  &
                       mat_reg(i,index)%ke_mat_d)
!
    ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Special case for Four by Four
!
         last = first + 3
         CALL v_mat_v_4(v(first:last,1:nj),                     &
                        v_scr(first:last,1:nj),                 &
                        mat_reg(i,index)%ke_mat_d)
    ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Special case for Five by Five
!
         last = first + 4
         CALL v_mat_v_5(v(first:last,1:nj),                     &
                        v_scr(first:last,1:nj),                 &
                        mat_reg(i,index)%ke_mat_d)
    ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Special case for Six by Six
!
         last = first + 5
         CALL v_mat_v_6(v(first:last,1:nj),                     &
                        v_scr(first:last,1:nj),                 &
                        mat_reg(i,index)%ke_mat_d)
    ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Special case for Seven by Seven
!
         last = first + 6
         CALL v_mat_v_7(v(first:last,1:nj),                     &
                        v_scr(first:last,1:nj),                 &
                        mat_reg(i,index)%ke_mat_d)
    ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Special case for Eight by Eight
!
         last = first + 7
         CALL v_mat_v_8(v(first:last,1:nj),                     &
                        v_scr(first:last,1:nj),                 &
                        mat_reg(i,index)%ke_mat_d)
    ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Special case for Nine by Nine
!
         last = first + 8
         CALL v_mat_v_9(v(first:last,1:nj),                     &
                        v_scr(first:last,1:nj),                 &
                        mat_reg(i,index)%ke_mat_d)
    ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Special case for Ten by Ten
!
         last = first + 9
         CALL v_mat_v_10(v(first:last,1:nj),                    &
                         v_scr(first:last,1:nj),                &
                         mat_reg(i,index)%ke_mat_d)
    END IF
    first = last
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE dvr_mat_mul_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!
!*deck dvr_mat_mul_2d_d
!***begin prologue     dvr_mat_mul_2d_d     
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
!***end prologue       dvr_mat_mul_2d_d
!
  SUBROUTINE dvr_mat_mul_2d_d(v,v_scr,ni,nj,nv,index)
  USE dvrprop_global,  ONLY   : mat_reg
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, nv, index
  REAL*8, DIMENSION(ni,nj,nv)              :: v
  REAL*8, DIMENSION(ni,nj,nv)              :: v_scr
  INTEGER                                  :: i, j
  INTEGER                                  :: first, last
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  first = 1
!
  DO i=1, num_reg(index)
     IF ( nfun_reg(i,index) > 10 ) then
!
!                        General Case
!
          last = first + nfun_reg(i,index) - 1
          DO j=1,nv
             CALL v_v_mat_gen(v(1:ni,first:last,j),             &
                              v_scr(1:ni,first:last,j),         &
                              mat_reg(i,index)%ke_mat_d) 
          END DO
!
     ELSE IF (nfun_reg(i,index) == 2) then
!
!                       Two by Two
! 
          last = first + 1
          DO j=1,nv
             CALL v_v_mat_2(v(1:ni,first:last,j),               &
                            v_scr(1:ni,first:last,j),           &
                            mat_reg(i,index)%ke_mat_d)
          END DO
     ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Three by Three
!
          last = first + 2
          DO j=1,nv
             CALL v_v_mat_3(v(1:ni,first:last,j),               &
                            v_scr(1:ni,first:last,j),           &
                            mat_reg(i,index)%ke_mat_d)

          END DO
     ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Four by Four
!
          last = first + 3
          DO j=1,nv
             CALL v_v_mat_4(v(1:ni,first:last,j),               &
                            v_scr(1:ni,first:last,j),           &
                            mat_reg(i,index)%ke_mat_d)
          END DO
     ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Five by Five
!
          last = first + 4
          DO j=1,nv
             CALL v_v_mat_5(v(1:ni,first:last,j),               &
                            v_scr(1:ni,first:last,j),           &
                            mat_reg(i,index)%ke_mat_d)
          END DO
     ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Six by Six
!
          last = first + 5
          DO j=1,nv
             CALL v_v_mat_6(v(1:ni,first:last,j),               &
                            v_scr(1:ni,first:last,j),           &
                            mat_reg(i,index)%ke_mat_d)
          END DO
     ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Seven by Seven
!
          last = first + 6
          DO j=1,nv
             CALL v_v_mat_7(v(1:ni,first:last,j),               &
                            v_scr(1:ni,first:last,j),           &
                            mat_reg(i,index)%ke_mat_d)
          END DO
     ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Eight by Eight
!
          last = first + 7
          DO j=1,nv
             CALL v_v_mat_8(v(1:ni,first:last,j),               &
                            v_scr(1:ni,first:last,j),           &
                            mat_reg(i,index)%ke_mat_d)
          END DO
     ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Nine by Nine
!
          last = first + 8
          DO j=1,nv
             CALL v_v_mat_9(v(1:ni,first:last,j),               &
                            v_scr(1:ni,first:last,j),           &
                            mat_reg(i,index)%ke_mat_d)
          END DO
     ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Ten by Ten
!
          last = first + 9
          DO j=1,nv
             CALL v_v_mat_10(v(1:ni,first:last,j),              &
                             v_scr(1:ni,first:last,j),          &
                             mat_reg(i,index)%ke_mat_d)
          END DO
     END IF
     first = last
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE dvr_mat_mul_2d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_diag_mul_1d_d
!***begin prologue     v_diag_mul_1d_d    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description         
!***                   
!***references
!***routines called    
!
!***end prologue       v_diag_mul_1d_d
!
  SUBROUTINE v_diag_mul_1d_d(v_in,v_out)
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1))             :: v_in
  REAL*8, DIMENSION(nphy(1))             :: v_out
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  v_out(:) = v_out(:) + grid(1)%v(:) * v_in(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END SUBROUTINE v_diag_mul_1d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_diag_mul_2d_d
!***begin prologue     v_diag_mul_2d_d    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description         
!***                   
!***references
!***routines called    
!
!***end prologue       v_diag_mul_2d_d
!
  SUBROUTINE v_diag_mul_2d_d(v_in,v_out)
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1))     :: v_in
  REAL*8, DIMENSION(nphy(2),nphy(1))     :: v_out
  INTEGER                                :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i = 1, nphy(2)
     v_out(i,:) = v_out(i,:) + grid(2)%v(i) * v_in(i,:)
  END DO
  DO i = 1, nphy(1)
     v_out(:,i) = v_out(:,i) + grid(1)%v(i) * v_in(:,i)
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END SUBROUTINE v_diag_mul_2d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_diag_mul_3d_d
!***begin prologue     v_diag_mul_3d_d    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description         
!***                   
!***references
!***routines called    
!
!***end prologue       v_diag_mul_3d_d
!
  SUBROUTINE v_diag_mul_3d_d(v_in,v_out)
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))     :: v_in
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))     :: v_out
  INTEGER                                        :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i = 1, nphy(3)
     v_out(i,:,:) = v_out(i,:,:) + grid(3)%v(i) * v_in(i,:,:)
  END DO
  DO i = 1, nphy(2)
     v_out(:,i,:) = v_out(:,i,:) + grid(2)%v(i) * v_in(:,i,:)
  END DO
  DO i = 1, nphy(1)
     v_out(:,:,i) = v_out(:,:,i) + grid(1)%v(i) * v_in(:,:,i)
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END SUBROUTINE v_diag_mul_3d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            END       MODULE dvr_matrix_vector_multiply_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
