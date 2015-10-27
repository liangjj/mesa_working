! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width_t=\linewidth_t}
! \title{V_tim}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck v_tim.f
!***begin prologue     v_tim
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           1-d hamiltonian
!***author             schneider, barry (nsf)
!***source
!***purpose            calculate one body time potential.
!***
!***
!***references
!***routines called
!***end prologue       v_tim
  SUBROUTINE v_tim
  USE prop_global
  USE Iterative_Global
  USE input_output
  IMPLICIT NONE
  INTEGER            :: first 
  REAL*8             :: fpkey
  CHARACTER (LEN=80) :: typ
  CHARACTER (LEN=80) :: chrkey
  LOGICAL            :: dollar, logkey, fretyp
  SAVE typ, fretyp
  DATA first / 0 /
  IF(first == 0) THEN
     IF(dollar(keywrd,card,cpass,inp)) THEN
        typ=chrkey(card,'potential','none',' ')
        fretyp=logkey(card,'convert_frequency',.false.,' ')
     END IF
     IF(units == 'atomic_units') THEN
        hbar = 1.d0
     END IF
     first = 1
  END IF
  IF(typ == 'none') THEN
     vt = 0.d0
  ELSE IF(typ == 't') THEN
     WRITE(iout,1)
     vt = t0
  ELSE IF(typ == 'cosine') THEN
     omega_t=fpkey(card,'electric_field_frequency',10.d0,' ')
     omega_t=omega_t/timau
     scale_t=fpkey(card,'scale_of_cosine_pulse',1.d0,' ')
     IF(fretyp) then
        omega_t = 2.d0 * pi * omega_t
     END IF
     IF(units == 'atomic_units') THEN
        omega_t = omega_t * timau
     END IF
     WRITE(iout,2) omega_t, scale_t
     vt = scale_t * COS(omega_t*t0)
  ELSE IF(typ == 'gaussian_pulse') THEN
     omega_t=fpkey(card,'electric_field_frequency',10.d0,' ')
     omega_t=omega_t/timau
     IF(fretyp) then
        omega_t = 2.d0 * pi * omega_t
     END IF
     IF(units == 'atomic_units') THEN
        omega_t = omega_t * timau
     END IF
     width_t=fpkey(card,'width_of_gaussian_pulse',0.d0,' ')
     shift_t=fpkey(card,'shift_of_gaussian_pulse',0.d0,' ')
     scale_t=fpkey(card,'scale_of_gaussian_pulse',1.d0,' ')
     WRITE(iout,3)
     WRITE(iout,4) scale_t, omega_t, width_t, shift_t
     vt = scale_t * COS( omega_t * t0 ) *  &
          EXP( - width_t * ( t0 - shift_t ) * ( t0 - shift_t ) )
  END IF
!  write(iout,*) vt
1    FORMAT(/,5X,'constructing factors for a time-dependent '  &
                 'perturbation = t')
2    FORMAT(/,5X,'constructing factors for a time-dependent '  &
                 'perturbation = A * cosine(omega*t)', &
            /,5X,'           omega = ',e15.8, &
            /,5X,'           A     = ',e15.8)
3    FORMAT(/,5X,'the form of the potential is:', &
            ///,10X,'v = scale * cos(omega*t) * ' &
                    'exp(-width*(t-shift)*(t-shift))')
4    FORMAT(/,5X,'constructing factors for a time-dependent '  &
                 'perturbation = A * cos(omega*t) * '  &
    '             exp(-W*(t-S)*(t-S))', &
            /,5X,'A        = ',e15.8, & 
            /,5X,'omega    = ',e15.8, &
            /,5X,'W        = ',e15.8, &
            /,5X,'S        = ',e15.8)
END SUBROUTINE v_tim




