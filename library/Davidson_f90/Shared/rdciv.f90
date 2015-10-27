!deck read_ci_vectors.f
!***begin prologue     read_ci_vectors
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           read ci vectors
!***author             schneider, barry (nsf)
!***source
!***
!***references
!***routines called
!***end prologue       read_ci_vectors
  SUBROUTINE read_ci_vectors(ci_vectors,n_root)
  IMPLICIT NONE
  REAL*8, DIMENSION(matrix_size,n_root)    :: ci_vectors
  INTEGER                                  :: n_root
  CHARACTER (LEN=4) :: itoc
  DO  i=1,nroot
      CALL iosys('read real "'//code//itoc(i)//'" from rwf', n,vec(1,i),0,' ')
!     write(iout,1) i, (vec(j,i),j=1,n)
END DO
RETURN
1    FORMAT(/,1X,'ci vector = ',i3,/,5E15.8)
END SUBROUTINE rdciv

