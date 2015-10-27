!deck add_absorb.f
!***begin prologue     add_absorb
!***date written       040706   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            addition of an absorbing potential to the 
!***                   kinetic energy matrix
!***references
!***routines called
!***end prologue       add_absorb
!
  SUBROUTINE add_absorb(pt,k_mat_d,v_add,k_mat_z,nr)
  INTEGER                                :: nr
  INTEGER                                :: i
  REAL*8, DIMENSION(nr,nr)               :: k_mat_d
  COMPLEX*16, DIMENSION(nr,nr)           :: k_mat_z
  COMPLEX*16                             :: eye=(0.d0,1.d0)
  REAL*8, DIMENSION(nr)                  :: pt
  COMPLEX*16, DIMENSION(nr)              :: v_add
!
!    Add an absorbing potential to the last element and store in a complex
!    kinetic energy matrix.
!    
  v_add = 1.d0
  DO i=1, nr
     DO j=1, nr
        k_mat_z = k_mat_d
     END DO
  END DO
  DO i=1,nr
     k_mat_z(i,i) = k_mat_z(i,i) - eye * v_add(i)
  END DO
END SUBROUTINE add_absorb
