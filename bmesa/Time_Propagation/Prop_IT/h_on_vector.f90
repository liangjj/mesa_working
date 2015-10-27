!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             MODULE h_on_vector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**begin prologue     h_on_vector
!**date written       040902   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            module for hamiltonian times vector.
!**description        this module contains all the routines necessary to
!**                   multiply a dvr Hamiltonian or a (3,5,7) point FD
!**                   Hamiltonian on a vector
!
!                     H = H(1,1) + H(2,2) + H(3,3)
!
!     Note that vectors are stored as V(nphy(3),nphy(2),nphy(1)).  Thus 
!     matrix vector multiplies in two and three dimensions must be done 
!     carefully for nphy(2) and nphy(3).  We have( repeated indices 
!     summed over), in 3D,
!    
!     V(3,2,1) =            H(1,1)   * V(3,2,1) 
!                                    + 
!                           H(2,2)   * V(3,2,1) 
!                                    + 
!                           H(3,3)   * V(3,2,1)   
!                                    =
!                           V(3,2,1) * H(1,1) 
!                                    + 
!                           V(3,2,1) * H(2,2) 
!                                    + 
!                           H(3,3) * V(3,2,1) 
!
! Thus the first mutiply may be done as a matrix V(3*2,1) on H(1,1), 
! the second as V(3,2,1) * H(2,2) with an outer loop over index 1 and 
! the third as a simple matrix vector multiply, H(3,3) * V(3,2*1).
!
!                             and in 2D,
!
!     V(2,1)   =            H(1,1) * V(2,1) 
!                                  + 
!                           H(2,2) * V(2,1)  
!                                  =
!                           V(2,1) * H(1,1) 
!                                  + 
!                           H(2,2) * V(2,1) 
!
!
!        
!**references
!**                   
!**routines called    
!**end prologue       h_on_vector
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!               The h_v and dvr_m_v routines perform the matrix
!               vector multiplication using a packed form of the matrix.
!               The finite_element_m_v routines use the element matrices
!               directly and have no need to pack the Hamiltonian.
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             CONTAINS
!deck h_v_d
!**begin prologue     h_v_d
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector in one dimension
!**
!**references
!**routines called h_m_v_1d_d, h_m_v_2d_d, h_m_v_3d_d
!**end prologue       h_v_d
  SUBROUTINE h_v_d(v_in,v_out,nv)
  USE h_m_v
  IMPLICIT NONE
  INTEGER                        :: nv, i
  REAL*8, DIMENSION(n3d,nv)      :: v_in          
  REAL*8, DIMENSION(n3d,nv)      :: v_out          
  v_out = 0.d0
!
  IF(spdim == 1) THEN
     CALL h_m_v_1d_d(v_in,v_out,nv)
  ELSE IF(spdim == 2) THEN
     CALL h_m_v_2d_d(v_in,v_out,nv)
  ELSE IF(spdim == 3) THEN
     CALL h_m_v_3d_d(v_in,v_out,nv)
  END IF  
!
!
! take care of the diagonal v
!
  DO i=1,n3d
     v_out(i,1:nv) = v_out(i,1:nv) + v_tot(i) * v_in(i,1:nv)
  END DO
  IF(log_prp(3)) THEN
     title='hamiltonian on vectors'
     CALL prntrm(title,v_out,n3d,nv,n3d,nv,iout)
  END IF
END SUBROUTINE h_v_d
!deck h_v_z
!**begin prologue     h_v_z
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector in one dimension
!**
!**references
!**routines called h_m_v_1d_z, h_m_v_2d_z, h_m_v_3d_z
!**end prologue       h_v_z
  SUBROUTINE h_v_z(v_in,v_out,nv)
  USE h_m_v
  IMPLICIT NONE
  INTEGER                                :: nv, i
  COMPLEX*16, DIMENSION(n3d,nv)          :: v_in          
  COMPLEX*16, DIMENSION(n3d,nv)          :: v_out          
  v_out = (0.d0,0.d0)
  IF(spdim == 1) THEN
     CALL h_m_v_1d_z(v_in,v_out,nv)
  ELSE IF(spdim == 2) THEN
     CALL h_m_v_2d_z(v_in,v_out,nv)
  ELSE IF(spdim == 3) THEN
     CALL h_m_v_3d_z(v_in,v_out,nv)
  END IF  
!
!
! take care of the diagonal v
!
  DO i=1,n3d
     v_out(i,1:nv) = v_out(i,1:nv) + v_tot(i) * v_in(i,1:nv)
  END DO
  IF(log_prp(3)) THEN
     title='hamiltonian on vectors'
     CALL prntcm(title,v_out,n3d,nv,n3d,nv,iout)
  END IF
END SUBROUTINE h_v_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             END MODULE h_on_vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
