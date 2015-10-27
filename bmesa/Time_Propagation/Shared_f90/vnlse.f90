!deck vnlse.f
!***begin prologue     vnlse
!***date written       041031   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           nonlinear schroediner equation
!***author             schneider, barry (nsf)
!***source
!***purpose            generate scaled parameters for nonlinear
!***                   schroedinger equation.  
!***references
!***routines called
!***end prologue       vnlse
!\begin{eqnarray}
!\end{eqnarray}
  SUBROUTINE vnlse
  USE potential
  USE dvr_global
  IMPLICIT NONE
  REAL*8                                 :: fpkey, tmp 
  INTEGER                                :: i, j, k, ii
  LOGICAL                                :: dollar
  IF( dollar( '$nlse',card,cpass,inp) )  THEN
      call fparr(card,'oscillator-frequencies-in-hertz',    &
                 omega_q,spdim,' ')
      natoms = fpkey(card,'number-of-atoms',1.d+06,' ')
      mass = fpkey(card,'mass',3.8192d-26,' ')
      scattering_length = fpkey(card,'scattering-length',   &
                                2.79d-09,' ')
!
!     Order them smallest to largest
!
      IF ( spdim > 1 ) THEN
           DO  ii=2,spdim
               i=ii-1
               k=i
               tmp=omega_q(i)
               DO  j=ii,spdim
                   IF(omega_q(j) < tmp) THEN
                      k=j
                      tmp=omega_q(j)
                   END IF
               END DO
               IF(k /= i) THEN
                  omega_q(k) = omega_q(i)
                  omega_q(i) = tmp
               END IF
           END DO
      END IF
  END IF
!
!     convert from Hertz to radians/sec
!
  omega_q = two * pi * omega_q
  DO i=2,spdim
     omega_q(i) = omega_q(i)/omega_q(1)
  END DO
!
!     compute fundamental scales
!
  d_0 = sqrt( hbar/(mass * omega_q(1)))
  e_0 = hbar * omega_q(1)
  t_0 = one/omega_q(1)
!
!     compute the coefficient of the nonlinear term
!
  nl_coef = four * pi * natoms * scattering_length/d_0

  write(iout,1) mass
  write(iout,2) (omega_q(i),i=1,spdim)
  write(iout,3) d_0, e_0, t_0, nl_coef
  write(inp,4) 
  write(inp,5) omega_q(1)
  write(inp,6)
  IF (spdim>1) THEN
      write(inp,7)
      write(inp,5) omega_q(2)
      write(inp,6)
  END IF
  IF (spdim>2) THEN
      write(inp,8)
      write(inp,5) omega_q(3)
      write(inp,6)
  END IF
 1    format(/,1x,'oscillator mass in kg      = ',e15.8)
 2    format(/,1x,'reduced oscillator-frequencies = ',3(e15.8,1x))
 3    format(/,1x,'length scale = ',e15.8,1x,'energy scale = ',e15.8,   &
             /,1x,'time scale   = ',e15.8,1x,'coupling constant = ',e15.8)
 4    format('$v_reg_1(x)')
 5    format(' potential=harmonic-oscillator omega=',d15.8)
 6    format('$end')
 7    format('$v_reg_1(y)')
 8    format('$v_reg_1(z)')
END SUBROUTINE vnlse
