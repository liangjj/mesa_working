!deck mwell.f
!***begin prologue     mwell
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            potential well
!***
!***references
!***routines called
!***end prologue       mwell
  SUBROUTINE mwell(v,pt,rwell,swell,n,nwell,prn)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n, nwell
  REAL*8, DIMENSION(n)                   :: v, pt
  REAL*8, DIMENSION(nwell+1)             :: rwell
  REAL*8, DIMENSION(nwell)               :: swell
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  INTEGER                                :: i, j
  DO i=1,nwell
     write(iout,1) rwell(i), rwell(i+1), swell(i)
  END DO
  DO  i=1,n
      DO j=1,nwell
         IF(pt(i) >= rwell(j).and.pt(i) <= rwell(j+1)) then
            v(i) = swell(j)
         ENDIF
      END DO
  END DO
  IF(prn) THEN
     title='potential'
     CALL prntrm(title,v,n,1,n,1,iout)
  END IF
1 FORMAT(/,5x,'value of potential between r         = ',e15.8, &
         /,5x,'           and                         '        &
         /,5x,'                           r         = ',e15.8, &
         /,5x,'                           potential = ',e15.8)  
END SUBROUTINE mwell



