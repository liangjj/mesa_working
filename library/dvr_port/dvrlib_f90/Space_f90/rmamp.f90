!deck rmamp
!***begin prologue     rmamp
!***date written       000623   (yymmdd)
!***revision date               (yymmdd)
!***keywords           r-matrix amplitudes for DVR basis
!***author             schneider, barry (nsf)
!***source
!***purpose            construct r-matrix amplitudes
!***description
!***references
!***routines called
!***end prologue       rmamp
  SUBROUTINE rmamp(eigvec,srf_0,srf,nphy)
  USE dvr_global,   ONLY   : output
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: nphy, i
  REAL*8, DIMENSION(nphy,nphy)           :: eigvec
  REAL*8, DIMENSION(2)                   :: srf_0
  REAL*8, DIMENSION(nphy,2)              :: srf
  CHARACTER (LEN=80) :: title
  srf(:,1) = eigvec(1,:)*srf_0(1)
  srf(:,2) = eigvec(nphy,:)*srf_0(2)
  IF(prn(11)) THEN
     title='r-matrix amplitudes'
     CALL prntrm(title,srf,nphy,2,nphy,2,output)
  END IF
END SUBROUTINE rmamp
