!deck wddtot.f
!***begin prologue     wddtot
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           eigenvalues, eigenvectors
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***                   
!***                   
!***                   
!***
!***description
!***references
!***routines called
!***end prologue       wddtot
  FUNCTION WDDOT(n,v_a,v_b,wt,pt)
  USE accuracy
  IMPLICIT NONE
  INTEGER                                :: n
  REAL(idp), DIMENSION(:)                :: v_a
  REAL(idp), DIMENSION(:)                :: v_b
  REAL(idp), DIMENSION(:), OPTIONAL      :: pt
  REAL(idp), DIMENSION(:)                :: wt
  REAL(idp)                              :: WDDOT
  INTEGER                                :: i
  WDDOT = 0.D0
  IF (present(pt)) THEN
      DO i=1,n
         WDDOT = WDDOT + v_a(i) * wt(i) * pt(i) * v_b(i)
      END DO 
  ELSE
      DO i=1,n
         WDDOT = WDDOT + v_a(i) * wt(i) * v_b(i)
      END DO 
  END IF
END FUNCTION WDDOT
