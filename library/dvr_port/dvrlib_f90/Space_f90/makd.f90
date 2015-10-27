!deck makd.f
!***begin prologue     makd
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            matrix elements between functions all in
!***                   region i.
!***references
!***routines called
!***end prologue       makd
!
  SUBROUTINE makd(hmat,kemat,keadd,pmat,pemat,peadd,norm,bridge,nr, &
                  nfun,nglobal,prn)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: nr, nfun, nglobal
  INTEGER                                :: i, j
  REAL*8, DIMENSION(nglobal,*)           :: hmat, pmat
  REAL*8, DIMENSION(nr,nr)               :: kemat, pemat
  REAL*8                                 :: keadd, peadd
  REAL*8, DIMENSION(nfun)                :: norm
  LOGICAL                                :: bridge, prn
  CHARACTER (LEN=80)                     :: title
!
  DO  i=1,nfun
      DO  j=1,i
          hmat(i,j) = hmat(i,j)  + norm(i)*norm(j)*kemat(i,j)
          pmat(i,j) = pmat(i,j)  + norm(i)*norm(j)*pemat(i,j)
      END DO
  END DO

!     we need to correct the last element if that element involves
!     a bridge function.

  IF(bridge) THEN
     hmat(nfun,nfun) = hmat(nfun,nfun) + norm(nfun)*norm(nfun)*keadd
     pmat(nfun,nfun) = pmat(nfun,nfun) + norm(nfun)*norm(nfun)*peadd
  END IF
  DO  i=1,nfun
      DO  j=1,i
          hmat(j,i) = hmat(i,j)
          pmat(j,i) = pmat(i,J)
      END DO
  END DO
  IF(prn) THEN
     title='matrix elements in region i'
     CALL prntrm(title,hmat,nfun,nfun,nglobal,nglobal,output)
  END IF
END SUBROUTINE makd



