!*********************************************************************
                      MODULE Asymptotic_Module
                      USE Coulomb_Variables_Module
!
                      IMPLICIT NONE
  REAL*8                                   :: eta_fac_plus
  REAL*8                                   :: eta_fac_minus
  REAL*8                                   :: eta_sq
  REAL*8                                   :: num_fac
  INTEGER                                  :: i_1
  INTEGER                                  :: i_2
  INTEGER                                  :: i_fac
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        Contains
!=======================================================================
!=======================================================================
!deck positive_energy_long_range_coefficients
!***begin prologue     long_range_coefficients
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            expansion coefficients for coulomb functions
!***                   at large rho.
!***description
!***references         nbs mathematical handbook
!***routines called
!***end prologue       long_range_coefficients
  SUBROUTINE positive_energy_long_range_coefficients(a_i,b_i)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                     :: a_i(0 : asymptotic_size)
  REAL*8, DIMENSION(:)                     :: b_i(0 : asymptotic_size)
  INTEGER                                  :: i
!**********************************************************************c
!                coefficients for asymptotic expansion                 c
!                of regular and irregular positive energy              c
!                        coulomb functions                             c
!**********************************************************************c
  eta_sq = eta_in*eta_in
  num_fac = angular_momentum * (angular_momentum + 1) + eta_sq
  DO i = 0, asymptotic_size
    i_1 = i + i + 1
    i_2 = i_1 + 1
    i_fac = i * (i + 1)
    a_i(i) = i_1 * eta_in / i_2
    b_i(i) = (num_fac - i_fac) / i_2
  END DO
  IF (print_long_range_coefficients) THEN
      WRITE(iout,*)
      WRITE(iout,*) '          the a_i and b_i long range expansion '//      &
                    'coefficients for L = '//itoc(angular_momentum)
      WRITE(iout,*)
      WRITE(iout,*) '      i          a_i               b_i'
      DO i = 0, asymptotic_size
         WRITE(iout,1) i, a_i(i), b_i(i)
      END DO
  END IF
1 FORMAT(5X,i3,3X,e15.8,3X,e15.8)
  END SUBROUTINE positive_energy_long_range_coefficients
!=======================================================================
!=======================================================================
!deck negative_energy_long_range_coefficients
!***begin prologue     long_range_coefficients
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            expansion coefficients for coulomb functions
!***                   at large rho.
!***description
!***references         nbs mathematical handbook
!***routines called
!***end prologue       long_range_coefficients
  SUBROUTINE negative_energy_long_range_coefficients(e_0)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                     :: e_0(0 : asymptotic_size)
  INTEGER                                  :: i
!**********************************************************************c
!             coefficients for asymtotic expansion of irregular        c
!             negative energy coulomb functions.                       c
!**********************************************************************c
  e_0(0) = one
  eta_fac_plus  = eta_in + angular_momentum
  eta_fac_minus = eta_in - angular_momentum - one
  DO i = 1 , asymptotic_size
     e_0(i) = ( i + eta_fac_plus ) * ( i + eta_fac_minus ) * e_0(i-1) / i
  END DO
  IF (print_long_range_coefficients) THEN
     WRITE(iout,*)
     WRITE(iout,*) '          the e_0 long range expansion '//      &
                    'coefficients for L = '//itoc(angular_momentum)
     WRITE(iout,*)
     WRITE(iout,*) '      i          e_0'
     DO i = 0, asymptotic_size
        WRITE(iout,1) i, e_0(i)
     END DO
  END IF
1 FORMAT(5X,i3,3X,e15.8)
  END SUBROUTINE negative_energy_long_range_coefficients
END MODULE Asymptotic_Module
