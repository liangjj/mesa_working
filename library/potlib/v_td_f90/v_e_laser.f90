!deck v_e_laser.f
!***begin prologue     v_e_laser
!***date written       050508   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            laser-electron interaction matrix
!***
!***description        calculate the time-dependent potential
!***                   for an electron interacting with a 
!***                   monchromatic laser.
!***
!***references

!***routines called
!***end prologue       v_e_laser
  SUBROUTINE v_e_laser(v,pt,e_field,t,a_0,omega,n_r,n)
  USE input_output
  IMPLICIT NONE
  REAL*8                                 :: pi=3.141592653589793238462643D+00
  REAL*8                                 :: c=137.03599976D+00
  INTEGER                                :: n_r, n
  REAL*8, DIMENSION(n)                   :: pt, v
  REAL*8                                 :: t, a_0
  REAL*8                                 :: omega
  REAL*8                                 :: switch, fac, s_fac, d_s_fac
  REAL*8                                 :: a_field, d_a_field, e_field
  switch=2.d0*pi*n_r/omega
  IF(t.ge.switch) THEN
     a_field =a_0
     e_field = - ( omega * a_field * cos(omega * t ) )/c
  ELSE
     fac = omega / ( 4.d0 * n_r )
     s_fac = sin( fac * t )
     d_s_fac = fac * cos ( fac * t )    
     a_field = a_0 * s_fac * s_fac
     d_a_field = 2.d0 * a_0 * s_fac * d_s_fac
     e_field = - ( d_a_field * sin (omega * t)             &
                           +                               &
                 omega * a_field * cos (omega * t) ) / c
  END IF
  v = e_field * pt
1    FORMAT(/,5X,'interaction potential = ',a32)
END SUBROUTINE v_e_laser














