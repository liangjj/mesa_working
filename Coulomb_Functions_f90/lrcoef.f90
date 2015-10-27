!deck long_range_coefficients

!***begin prologue     long_range_coefficients
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            expansion coefficients needed to construct
!***                   asymptotic values for regular or irregular
!***                   coulomb functions at large rho.
!***
!***references         nbs handbook

!***routines called

!***end prologue       long_range_coefficients
  SUBROUTINE long_range_coefficients(a_i,b_i,e_0)
  IMPLICIT NONE
  REAL*8, DIMENSION(0 : asy                      :: ak(0:n)
  REAL*8                      :: bk(0:n)
  REAL*8                      :: e0(0:n)
  REAL*8                      :: eta_sq, numfac
  REAL*8                      :: num_fac


COMMON/io/inp,iout

IF (energy == 'positive') THEN
!**********************************************************************c
!                coefficients for asymptotic expansion                 c
!                of regular and irregular positive energy              c
!                        coulomb functions                             c
!**********************************************************************c
  etasq=eta*eta
  numfac=l*(l+1)+etasq
  DO  k=0,n
    k1=k+k+1
    k2=k1+1
    kfac=k*(k+1)
    ak(k)=k1*eta/k2
    bk(k)=(numfac-kfac)/k2
  END DO
  IF (level >= 2) THEN
    WRITE(iout,*)
    WRITE(iout,*) '          the a and b long range '//  &
        'expansion coefficients'
    WRITE(iout,*)
    WRITE(iout,*) '      i          a                 b'
    DO  i=0,n
      WRITE(iout,1) i,ak(i), bk(i)
    END DO
  END IF
ELSE IF (energy == 'negative') THEN
!**********************************************************************c
!             coefficients for asymtotic expansion of irregular        c
!             negative energy coulomb functions.                       c
!**********************************************************************c
  e0(0)=1.d0
  DO  i=1,n
    e0(i)=(i+eta+l)*(i+eta-l-1)*e0(i-1)/i
  END DO
  IF (level >= 2) THEN
    WRITE(iout,*)
    WRITE(iout,*) '          the e long range expansion '// 'coefficients'
    WRITE(iout,*)
    WRITE(iout,*) '      i          e'
    DO  i=0,n
      WRITE(iout,2) i,e0(i)
    END DO
  END IF
END IF
RETURN
1 FORMAT(5X,i3,3X,e15.8,3X,e15.8)
2 FORMAT(5X,i3,3X,e15.8)
END SUBROUTINE long_range_coefficients













NIST, Gaithersburg, MD    Anne.Chaka@nist.gov
