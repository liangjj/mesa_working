!deck ovmat.f
!***begin prologue     ovmat
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            lobatto overlap ovrix elements.
!***
!***references

!***routines called
!***end prologue       ovmat

  SUBROUTINE ovmat(ovr,far,fbr,wtr,nr,region)
  USE dvr_global ,  ONLY  : iout
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                              :: nr
  REAL*8, DIMENSION(nr,nr)             :: ovr, far, fbr
  REAL*8, DIMENSION(nr)                :: wtr
  INTEGER                              :: region
  CHARACTER (LEN=80)                   :: title
  CHARACTER (LEN=3)                    :: itoc
  INTEGER                              :: i
  CALL rzero(ovr,nr*nr)
  DO  i=1,nr
      ovr(i,i) = ovr(i,i) + far(i,i)*wtr(i)*fbr(i,i)
  END DO
  IF(prn(3)) THEN
     title='unnormalized overlap matrix for sector = '//itoc(region)
     CALL prntrm(title,ovr,nr,nr,nr,nr,iout)
  END IF
END SUBROUTINE ovmat
