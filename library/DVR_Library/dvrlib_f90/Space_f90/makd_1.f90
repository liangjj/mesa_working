!deck makd_1.f
!***begin prologue     makd_1
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            matrix elements between functions all in
!***                   region i.
!***references
!***routines called
!***end prologue       makd_1
!
  SUBROUTINE makd_1(hmat,kemat,keadd,norm,bridge,nr, &
                  nfun,nglobal,prn)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: nr, nfun, nglobal
  INTEGER                                :: i, j
  REAL*8, DIMENSION(nglobal,*)           :: hmat
  REAL*8, DIMENSION(nr,nr)               :: kemat
  REAL*8                                 :: keadd
  REAL*8, DIMENSION(nfun)                :: norm
  LOGICAL                                :: bridge, prn
  CHARACTER (LEN=80)                     :: title
!
  DO  i=1,nfun
      DO  j=1,i
          hmat(i,j) = hmat(i,j)  + norm(i)*norm(j)*kemat(i,j)
      END DO
  END DO

!     we need to correct the last element if that element involves
!     a bridge function.

  IF(bridge) THEN
     hmat(nfun,nfun) = hmat(nfun,nfun) + norm(nfun)*norm(nfun)*keadd
  END IF
  DO  i=1,nfun
      DO  j=1,i
          hmat(j,i) = hmat(i,j)
      END DO
  END DO
  IF(prn) THEN
     title='matrix elements in region i'
     CALL prntrm(title,hmat,nfun,nfun,nglobal,nglobal,iout)
  END IF
END SUBROUTINE makd_1



