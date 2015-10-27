!=======================================================================
  Program test
!=======================================================================
!  Get information from the user about the problem to be solved
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    REAL(KIND=8) :: z
    CHARACTER(LEN=3) :: term

    write(*, '(A)',ADVANCE='NO') 'Atomic Number: '
    read(*,*) z
    write(*, '(A)',ADVANCE='NO') 'Term         : '
    read(*,*) term
    print *, z,term
  END program test
