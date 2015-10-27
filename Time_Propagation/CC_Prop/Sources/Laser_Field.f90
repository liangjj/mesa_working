!deck Laser_Field
!***begin prologue     Laser_Field
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose             Calculate electric field at time = t, where            
!***                    E(t) = E_0 f(t) sin(omega*t)
!***                           E_0 is the peak field strength
!***references
!***routines called    Name             Location
!                      ----             --------
!                    
!                    
!
!***modules used       Name             Location
!                      ----             --------
!                   
!                   
!                   
!
!***end prologue       Laser_Field  
  Subroutine Laser_Field (t)  
!============================================================== !
! Calculate the time-dependent electric field at a given time.  ! 
! 't' is in atomic units and 'pulse_duration' is given in units ! 
! of optical cycles.                                            !  
!============================================================== !
! Input  :: t, peak_electric_field, pulse_duration
! Output :: electric_field
                 USE grid_global,                    ONLY: pi, two_pi
                 USE Global_Time_Propagation_Module, ONLY:            &
                                                     photon_energy,   &
                                                     type_pulse,      &
                                                     electric_field,  &
                                                     pulse_duration,  &
                                                     ramp_time,       &
                                                     peak_electric_field
  IMPLICIT NONE
  REAL*8                         :: t
  REAL*8                         :: tau
  REAL*8                         :: top
  REAL*8                         :: temp
!
!           Convert to optical cycles
!
  tau = photon_energy * t /two_pi
  IF ( tau < 0.d0) THEN
       Call lnkerr('error: time improperly defined')
       stop
  END IF       
  IF (type_pulse == 'sine') THEN
!                      Sine^2 Profile 
      temp = SIN(tau * pi / pulse_duration)
      electric_field = peak_electric_field * temp * temp * cos(two_pi * tau)
  ELSE IF( type_pulse == 'flat_top') THEN
! 
          electric_field = 0.d0
          IF ( tau > pulse_duration) THEN
               return
          ELSE
               top = pulse_duration - ramp_time - ramp_time
               IF (tau < ramp_time) THEN
                   temp = tau / ramp_time
               ELSE IF ( tau >= ramp_time .and. tau < (pulse_duration - ramp_time) ) THEN
                   temp = 1.0D+00
               ELSE IF ( tau >= (pulse_duration - ramp_time) ) THEN
                   temp = (pulse_duration - tau) / ramp_time
               END IF
          END IF
          electric_field = peak_electric_field * temp * SIN(two_pi * tau)
  END IF
  END SUBROUTINE Laser_Field
