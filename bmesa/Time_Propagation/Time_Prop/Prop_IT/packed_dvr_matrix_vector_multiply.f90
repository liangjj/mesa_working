!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      MODULE packed_dvr_matrix_vector_multiply
                         USE arnoldi_global
                         USE dvrprop_global
                         USE dvr_shared
                         USE dvr_global
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck dvr_m_v_1d_d
!***begin prologue     dvr_m_v_1d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Multiply a packed DVR hamiltonian on a vector
!***                   There are other routines which do this operation
!**                    avoiding packing and recognizing the intrisic
!**                    structure of the FEDVR Hamiltonian.  This module
!**                    was written much earlier and is based on the
!**                    typical quantum chemistry procedure of separating
!**                    out the diagonal and then packing the non-zero
!**                    elements and their indices into arrays.  The
!**                    2D array hibuf contains the indices and the array
!**                    hbuf the matrix element.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       dvr_m_v_1d_d
  SUBROUTINE dvr_m_v_1d_d(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                :: nv
  REAL*8, DIMENSION(nphy(1),nv)          :: v_in          
  REAL*8, DIMENSION(nphy(1),nv)          :: v_out          
  REAL*8                                 :: hij
  INTEGER                                :: i, j, nel 
!
!    First add in the diagonal contribution.
!
  DO i=1,nphy(1)
     v_out(i,1:nv) = v_out(i,1:nv)                       &
                                     +                   &
                     buf(1)%d(i) * v_in(i,1:nv)
  END DO
!
! Now loop over the non-zero, off-diagonal eleemnts.
!
  DO  nel=1,nonz(1)
      i=buf(1)%hibuf(1,nel)
      j=buf(1)%hibuf(2,nel)
      hij=buf(1)%hbuf(nel)
      v_out(i,1:nv) = v_out(i,1:nv) + hij * v_in(j,1:nv)
      v_out(j,1:nv) = v_out(j,1:nv) + hij * v_in(i,1:nv)
  END DO
END SUBROUTINE dvr_m_v_1d_d
!deck dvr_m_v_1d_z
!***begin prologue     dvr_m_v_1d_z
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       dvr_m_v_1d_z
  SUBROUTINE dvr_m_v_1d_z(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                :: nv
  COMPLEX*16, DIMENSION(nphy(1),nv)      :: v_in          
  COMPLEX*16, DIMENSION(nphy(1),nv)      :: v_out          
  REAL*8                                 :: hij
  INTEGER                                :: i, j, nel 
!
!    First add in the diagonal contribution.
!
  DO i=1,nphy(1)
     v_out(i,1:nv) = v_out(i,1:nv)                     &
                                    +                  &
                     buf(1)%d(i) * v_in(i,1:nv)
  END DO
!
! Now loop over the non-zero, off-diagonal eleemnts.
!
  DO  nel=1,nonz(1)
      i=buf(1)%hibuf(1,nel)
      j=buf(1)%hibuf(2,nel)
      hij=buf(1)%hbuf(nel)
      v_out(i,1:nv) = v_out(i,1:nv) + hij * v_in(j,1:nv)
      v_out(j,1:nv) = v_out(j,1:nv) + hij * v_in(i,1:nv)
  END DO
END SUBROUTINE dvr_m_v_1d_z
!deck dvr_m_v_2d_d.f
!***begin prologue     dvr_m_v_2d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       dvr_m_v_2d_d
  SUBROUTINE dvr_m_v_2d_d(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                     :: nv
  REAL*8, DIMENSION(nphy(2),nphy(1),nv)       :: v_in, v_out
  REAL*8                                      :: hij
  INTEGER                                     :: i, nel, j
!
!                      Perform the multiply in nphy(2)
!
  DO  i=1,nphy(2)
      v_out(i,1:nphy(1),1:nv) = v_out(i,1:nphy(1),1:nv)              &
                                         +                           &
                                buf(2)%d(i) * v_in(i,1:nphy(1),1:nv)
  END DO
  DO  nel=1,nonz(2)
      i=buf(2)%hibuf(1,nel)
      j=buf(2)%hibuf(2,nel)
      hij=buf(2)%hbuf(nel)
      v_out(i,1:nphy(1),1:nv) = v_out(i,1:nphy(1),1:nv)              &
                                         +                           &
                               hij * v_in(j,1:nphy(1),1:nv)
      v_out(j,1:nphy(1),1:nv) = v_out(j,1:nphy(1),1:nv)              &
                                         +                           &
                                hij * v_in(i,1:nphy(1),1:nv)
  END DO
!
!                      Perform the multiply in nphy(1)
!
  DO  i=1,nphy(1)
      v_out(1:nphy(2),i,1:nv) = v_out(1:nphy(2),i,1:nv)              &
                                         +                           &
                                buf(1)%d(i) * v_in(1:nphy(2),i,1:nv)
  END DO
  DO  nel=1,nonz(1)
      i=buf(1)%hibuf(1,nel)
      j=buf(1)%hibuf(2,nel)
      hij=buf(1)%hbuf(nel)
      v_out(1:nphy(2),i,1:nv) = v_out(1:nphy(2),i,1:nv)              &
                                         +                           &
                               hij * v_in(1:nphy(2),j,1:nv)
      v_out(1:nphy(2),j,1:nv) = v_out(1:nphy(2),j,1:nv)              &
                                         +                           &
                                hij * v_in(1:nphy(2),i,1:nv)
  END DO
END SUBROUTINE dvr_m_v_2d_d
!deck dvr_m_v_2d_z.f
!***begin prologue     dvr_m_v_2d_z
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       dvr_m_v_2d_z
  SUBROUTINE dvr_m_v_2d_z(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                     :: nv
  COMPLEX*16, DIMENSION(nphy(2),nphy(1),nv)   :: v_in
  COMPLEX*16, DIMENSION(nphy(2),nphy(1),nv)   :: v_out
  REAL*8                                      :: hij
  INTEGER                                     :: i, nel, j
  DO  i=1,nphy(2)
      v_out(i,1:nphy(1),1:nv) = v_out(i,1:nphy(1),1:nv)                 &
                                         +                              &
                                buf(2)%d(i) * v_in(i,1:nphy(1),1:nv)
  END DO
  DO  nel=1,nonz(2)
      i=buf(2)%hibuf(1,nel)
      j=buf(2)%hibuf(2,nel)
      hij=buf(2)%hbuf(nel)
      v_out(i,1:nphy(1),1:nv) = v_out(i,1:nphy(1),1:nv)                  &
                                         +                               &
                                hij * v_in(j,1:nphy(1),1:nv)
      v_out(j,1:nphy(1),1:nv) = v_out(j,1:nphy(1),1:nv)                  &
                                         +                               &
                                hij * v_in(i,1:nphy(1),1:nv)
  END DO
  DO  i=1,nphy(1)
      v_out(1:nphy(2),i,1:nv) = v_out(1:nphy(2),i,1:nv)                  &
                                         +                               &
                                buf(1)%d(i) * v_in(1:nphy(2),i,1:nv)
  END DO
  DO  nel=1,nonz(1)
      i=buf(1)%hibuf(1,nel)
      j=buf(1)%hibuf(2,nel)
      hij=buf(1)%hbuf(nel)
      v_out(1:nphy(2),i,1:nv) = v_out(1:nphy(2),i,1:nv)                  &
                                               +                         &
                                hij * v_in(1:nphy(2),j,1:nv)
      v_out(1:nphy(2),j,1:nv) = v_out(1:nphy(2),j,1:nv)                  &
                                               +                         &
                                hij * v_in(1:nphy(2),i,1:nv)
  END DO
END SUBROUTINE dvr_m_v_2d_z
!deck dvr_m_v_3d_d.f
!***begin prologue     dvr_m_v_3d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       dvr_m_v_3d_d
  SUBROUTINE dvr_m_v_3d_d(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                             :: nv 
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),nv)       :: v_in, v_out
  REAL*8                                              :: hij
  INTEGER                                             :: i, nel, j
!
!                Perform the multiply in nphy(3)
!
  DO i=1,nphy(3)
     v_out(i,1:nphy(2),1:nphy(1),1:nv)                     & 
                         =                                 &
     v_out(i,1:nphy(2),1:nphy(1),1:nv)                     &
                         +                                 &
     buf(3)%d(i) * v_in(i,1:nphy(2),1:nphy(1),1:nv)     
  END DO
  DO nel=1,nonz(3)
      i=buf(3)%hibuf(1,nel)
      j=buf(3)%hibuf(2,nel)
      hij=buf(3)%hbuf(nel)
      v_out(i,1:nphy(2),1:nphy(1),1:nv)                    &
                        =                                  &
      v_out(i,1:nphy(2),1:nphy(1),1:nv)                    &
                        +                                  &
           hij * v_in(j,1:nphy(2),1:nphy(1),1:nv)
      v_out(j,1:nphy(2),1:nphy(1),1:nv)                    &
                        =                                  &
      v_out(j,1:nphy(2),1:nphy(1),1:nv)                    &
                        +                                  &
      hij * v_in(i,1:nphy(2),1:nphy(1),1:nv)
  END DO
!
!                Perform the multiply in nphy(2)
!
  DO  i=1,nphy(2)
      v_out(1:nphy(3),i,1:nphy(1),1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),i,1:nphy(1),1:nv)                    &
                        +                                  &
      buf(2)%d(i) * v_in(1:nphy(3),i,1:nphy(1),1:nv) 
  END DO
  DO  nel=1,nonz(2)
      i=buf(2)%hibuf(1,nel)
      j=buf(2)%hibuf(2,nel)
      hij=buf(2)%hbuf(nel)
      v_out(1:nphy(3),i,1:nphy(2),1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),i,1:nphy(1),1:nv)                    &
                        +                                  & 
              hij * v_in(1:nphy(3),j,1:nphy(1),1:nv)
      v_out(1:nphy(3),j,1:nphy(1),1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),j,1:nphy(1),1:nv)                    &
                        +                                  &
      hij * v_in(1:nphy(3),i,1:nphy(1),1:nv)
  END DO
!
!                Perform the multiply in nphy(1)
!
  DO  i=1,nphy(1)
      v_out(1:nphy(3),1:nphy(2),i,1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),1:nphy(2),i,1:nv)                    &
                        +                                  &
      buf(1)%d(i) * v_in(1:nphy(3),1:nphy(2),i,1:nv)
  END DO
  DO  nel=1,nonz(1)
      i=buf(1)%hibuf(1,nel)
      j=buf(1)%hibuf(2,nel)
      hij=buf(1)%hbuf(nel)
      v_out(1:nphy(3),1:nphy(2),i,1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),1:nphy(2),i,1:nv)                    &
                        +                                  &
             hij * v_in(1:nphy(3),1:nphy(2),j,1:nv)        
      v_out(1:nphy(3),1:nphy(2),j,1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),1:nphy(2),j,1:nv)                    & 
                        +                                  &
      hij * v_in(1:nphy(3),1:nphy(2),i,1:nv)
  END DO
END SUBROUTINE dvr_m_v_3d_d 
!deck dvr_m_v_3d_z.f
!***begin prologue     dvr_m_v_3d_z
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       dvr_m_v_3d_z
  SUBROUTINE dvr_m_v_3d_z(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                              :: nv
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1),nv)    :: v_in, v_out
  REAL*8                                               :: hij
  INTEGER                                              :: i, nel, j
  DO i=1,nphy(3)
     v_out(i,1:nphy(2),1:nphy(1),1:nv)                     & 
                         =                                 &
     v_out(i,1:nphy(2),1:nphy(1),1:nv)                     &
                         +                                 &
     buf(3)%d(i) * v_in(i,1:nphy(2),1:nphy(1),1:nv)     
  END DO
  DO nel=1,nonz(3)
      i=buf(3)%hibuf(1,nel)
      j=buf(3)%hibuf(2,nel)
      hij=buf(3)%hbuf(nel)
      v_out(i,1:nphy(2),1:nphy(1),1:nv)                    &
                        =                                  &
      v_out(i,1:nphy(2),1:nphy(1),1:nv)                    &
                        +                                  &
           hij * v_in(j,1:nphy(2),1:nphy(1),1:nv)
      v_out(j,1:nphy(2),1:nphy(1),1:nv)                    &
                        =                                  &
      v_out(j,1:nphy(2),1:nphy(1),1:nv)                    &
                        +                                  &
      hij * v_in(i,1:nphy(2),1:nphy(1),1:nv)
  END DO
  DO  i=1,nphy(2)
      v_out(1:nphy(3),i,1:nphy(1),1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),i,1:nphy(1),1:nv)                    &
                        +                                  &
      buf(2)%d(i) * v_in(1:nphy(3),i,1:nphy(1),1:nv) 
  END DO
  DO  nel=1,nonz(2)
      i=buf(2)%hibuf(1,nel)
      j=buf(2)%hibuf(2,nel)
      hij=buf(2)%hbuf(nel)
      v_out(1:nphy(3),i,1:nphy(2),1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),i,1:nphy(1),1:nv)                    &
                        +                                  & 
      hij * v_in(1:nphy(3),j,1:nphy(1),1:nv)
      v_out(1:nphy(3),j,1:nphy(1),1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),j,1:nphy(1),1:nv)                    &
                        +                                  &
      hij * v_in(1:nphy(3),i,1:nphy(1),1:nv)
  END DO
  DO  i=1,nphy(1)
      v_out(1:nphy(3),1:nphy(2),i,1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),1:nphy(2),i,1:nv)                    &
                        +                                  &
      buf(1)%d(i) * v_in(1:nphy(3),1:nphy(2),i,1:nv)
  END DO
  DO  nel=1,nonz(1)
      i=buf(1)%hibuf(1,nel)
      j=buf(1)%hibuf(2,nel)
      hij=buf(1)%hbuf(nel)
      v_out(1:nphy(3),1:nphy(2),i,1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),1:nphy(2),i,1:nv)                    &
                        +                                  &
             hij * v_in(1:nphy(3),1:nphy(2),j,1:nv)        
      v_out(1:nphy(3),1:nphy(2),j,1:nv)                    &
                        =                                  &
      v_out(1:nphy(3),1:nphy(2),j,1:nv)                    & 
                        +                                  &
      hij * v_in(1:nphy(3),1:nphy(2),i,1:nv)
  END DO
END SUBROUTINE dvr_m_v_3d_z 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            END       MODULE packed_dvr_matrix_vector_multiply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
