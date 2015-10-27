! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{V_couple}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck v_couple.f
!**begin prologue     v_couple
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            time-dependent potential
!**
!**description        input data for coupling potential, including the
!**                   mean field interaction for a contact potential.
!**references
!**routines called
!**end prologue       v_couple
  SUBROUTINE v_couple
  USE arnoldi_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                      :: natom, intkey
  REAL*8                       :: fpkey, sctlen
  CHARACTER (LEN=80)           :: chrkey
  CHARACTER (LEN=2)            :: atom
  LOGICAL                      :: dollar, fretyp, logkey
  IF(dollar('$v_couple',card,cpass,inp) ) THEN
     vtyp(1)=chrkey(card,'v-space-time','none',' ')
     vtyp(2)=chrkey(card,'v-non-linear','none',' ')
     IF(vtyp(2) /= 'none') THEN
        atom=chrkey(card,'atom','cs',' ')
        sctlen=fpkey(card,'scattering-length',52.d0,' ')
        natom=intkey(card,'number-of-atoms',0,' ')
     END IF
     WRITE(iout,1) vtyp
     WRITE(iout,2) units
     fretyp=logkey(card,'convert-frequency',.false.,' ')
     scale=fpkey(card,'electric-field-strength',1.d0,' ')
     scale=scale*efieldau
     omega=fpkey(card,'electric-field-frequency',1.d0,' ')
     omega=omega/timau
     IF(vtyp(1) == 'cosine-gaussian-pulse'.OR.  &
        vtyp(1) == 'sine-gaussian-pulse') THEN
        width=fpkey(card,'width-of-gaussian-pulse',0.d0,' ')
        shift=fpkey(card,'shift-of-gaussian-pulse',0.d0,' ')
     END IF
     IF(fretyp) then
        omega=2.d0*pi*omega
     END IF
     IF(units == 'atomic-units') THEN
        omega=omega*timau
        scale=scale/efieldau
     END IF
     IF(vtyp(1) == 'cosine-pulse') THEN
        WRITE(iout,3) omega, scale
     ELSE IF(vtyp(1) == 'sine-pulse') THEN
        WRITE(iout,4) omega, scale
     ELSE IF(vtyp(1) == 't') THEN
        WRITE(iout,5)
     ELSE IF(vtyp(1) == 'cosine-dipole-field') THEN
        WRITE(iout,6) omega, scale
     ELSE IF(vtyp(1) == 'sine-dipole-field') THEN
        WRITE(iout,7) omega, scale
     ELSE IF(vtyp(1) == 'cosine-gaussian-pulse') THEN
        WRITE(iout,8) scale, omega, width, shift
     ELSE IF(vtyp(1) == 'sine-gaussian-pulse') THEN
        WRITE(iout,9) scale, omega, width, shift
     ELSE IF(vtyp(1) == 'interferometer') THEN
        n_scale=intkey(card,'n_i',10,' ')
        e_c=fpkey(card,'e_c',.001,' ')
        tau_c=fpkey(card,'tau',5.d-03,' ')
        k_0=fpkey(card,'k_0',1.0,' ')
        write(iout,12) n_scale, e_c, tau_c, k_0 
     END IF
     IF(natom /= 0) THEN
        WRITE(iout,11) atom, natom, sctlen
        IF(atom == 'cs') THEN
           mass=2.2D-25
        ELSE IF(atom == 'na') THEN
           mass=3.8176D-26
        END IF
        IF(units == 'atomic-units') THEN
           hbar=1.d0
           mass=mass/massau
        END IF
        gamma = natom*natom*4.d0*pi*hbar*hbar * sctlen/(mass*mass)
     END IF
  ELSE
     write(iout,13)
     stop
  END IF
1    FORMAT(/,1X,'space and time    = ',a32, &
            /,1X, 'non-linear        = ',a32)
2    FORMAT(/,5X,'units = ',a16)
3    FORMAT(/,5X,'perturbation = E0 * cosine(omega * t)',  &
            /,5X,'       omega = ',e15.8, &
            /,5X,'       E0    = ',e15.8)
4    FORMAT(/,5X,'perturbation = E0 * sine(omega * t)',  &
            /,5X,'       omega = ',e15.8, &
            /,5X,'       E0    = ',e15.8)
5    FORMAT(/,5X,'perturbation =  t')
6    FORMAT(/,5X,'perturbation = E0 * x * cosine(omega * t)',  &
            /,5X,'       omega = ',e15.8, &
            /,5X,'       E0    = ',e15.8)
7    FORMAT(/,5X,'perturbation = E0 * x * sine(omega * t)',  &
            /,5X,'       omega = ',e15.8, &
            /,5X,'       E0    = ',e15.8)
8    FORMAT(/,5X,'perturbation = A * x * cos(omega * t) * ',  &
            /,5X,'               exp(- W * (t-S) * (t-S))',  &
            /,5X,'    A        = ',e15.8, &
            /,5X,'    omega    = ',e15.8, &
            /,5X,'    W        = ',e15.8, &
            /,5X,'    S        = ',e15.8)
9    FORMAT(/,5X,'perturbation = A * x * sin(omega * t) * ',  &
            /,5X,'               exp(- W * (t-S) * (t-S))',  &
            /,5X,'    A        = ',e15.8, &
            /,5X,'    omega    = ',e15.8, &
            /,5X,'    W        = ',e15.8, &
            /,5X,'    S        = ',e15.8)
11   FORMAT(/,15X,'atomic parameters',/,           &
            /,5X, 'atom                    = ',a2, &
            /,5X, 'number of atoms in trap = ',i8, &
            /,5X, 'scattering length       = ',e15.8)
12   FORMAT(/,5X,'perturbation = [N * k(t) / E_c ] * '       &
                 '[ 1 - cos( x ) ]'                          &
            /,5X,'    k(t)     = k(0) * exp ( - t / tau )'   &
            /,5X,'    N        = ',i6,                       &
            /,5X,'    E_c      = ',e15.8,                    &
            /,5X,'    tau      = ',e15.8,                    &
            /,5X,'    k(0)     = ',e15.8)
13   FORMAT(/,1x,'potential keyword error')
END SUBROUTINE v_couple


