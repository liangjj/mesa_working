!deck add_absorb.f
!***begin prologue     add_absorb
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            addition of an absorbing potential to kinetic energy
!***references
!***routines called
!***end prologue       add_absorb
!
  SUBROUTINE add_absorb(k_mat_d,k_mat_z,nr)
  USE io
  INTEGER                                :: nr
  INTEGER                                :: i
  REAL*8, DIMENSION(nr,nr)               :: k_mat_d
  COMPLEX*16, DIMENSION(nr,nr)           :: k_mat_z
  COMPLEX*16                             :: eye=(0.d0,1.d0)
  REAL*8                                 :: v_add=1.d-01
!
!    Add an absorbing potential to the last element and store in a complex
!    kinetic energy matrix.
!    
  DO i=1, nr
     DO j=1, nr
        k_mat_z = k_mat_d
     END DO
  END DO
  DO i=1,nr
     k_mat_z(i,i) = k_mat_z(i,i) - eye * v_add
  END DO
END SUBROUTINE add_absorb
