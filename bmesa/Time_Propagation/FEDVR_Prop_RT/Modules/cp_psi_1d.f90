!*deck cp_psi_1d
!**begin prologue     cp_psi_1d
!**date written       040707   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate the spatial part of the 
!                     initial wavepacket.
!**references
!**                   
!**routines called    c_vect, rad_packet, gauss_packet, sppose
!**end prologue       cp_psi_1d
!
  SUBROUTINE cp_psi_1d(wave_function,scratch_vector,t,t_calc)
  USE prop_global
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  USE normalize
  USE moment
  USE plot_wavefunction
  USE gaussian_wave_packet
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),2)           :: wave_function
  REAL*8, DIMENSION(nphy(1),2)           :: scratch_vector
  INTEGER                                :: t, i 
  REAL*8                                 :: norm, nfac
  REAL*8, DIMENSION(2)                   :: hfac
  LOGICAL                                :: dollar, logkey
  REAL*8                                 :: ddot
  REAL*8                                 :: t_calc
  CHARACTER (LEN=2)                      :: itoc
  CHARACTER (LEN=80)                     :: chrkey
  INTEGER                                :: intkey, state
  SAVE norm
  hfac(1)=1.d0/hbar
  hfac(2)=.5d0*hfac(1)
  IF(t == 1) THEN
     IF( dollar('$initial-state',card,title,inp) ) THEN
         i0stat=chrkey(card,'driver','unperturbed-state-vector',' ')
         prnton=logkey(card,'print=on',.false.,' ')
         WRITE(iout,1) t0, i0stat
     END IF
  
!        The allowable right hand sides are; one, a state vector of the
!        unperturbed Hamiltonian, a specified superposition of 
!        states of the unperturbed Hamiltonian or a Gaussian 
!        pulse with a velocity component.  The routine could easily be 
!        modified to do more.
     IF(i0stat == 'one') THEN
        wave_function(:,1) =    cos (t0*t0*hfac(2))
        wave_function(:,2) =  - sin (t0*t0*hfac(2))
     ELSE IF(i0stat.eq.'unperturbed-state-vector') then
        state=intkey(card,'initial-state',0,' ')
        energy=grid(1)%eigv_0(state)
        wave_function(:,1) = grid(1)%eigvec_0(:,state)   
!
!       Use scratch_vector as a temporary
!          
        scratch_vector(:,1) =   wave_function(:,1)                      &
                                    *                                   &
                                cos(energy*t0*hfac(1))   
        scratch_vector(:,2) = - wave_function(:,1)                      &
                                    *                                   &
                                sin(energy*t0*hfac(1)) 
!
!       Now store it
! 
        wave_function = scratch_vector
     ELSE IF(i0stat == 'perturbed-state-vector') then
        state=intkey(card,'initial-state',0,' ')
        energy=grid(1)%eigv(state+1)
        wave_function(:,1) =   cos(energy*t0*hfac(1))                   &
                                         *                              &
                               grid(1)%eigvec(:,state+1)
        wave_function(:,2) = - sin(energy*t0*hfac(1))                   &
                                         *                              &
                               grid(1)%eigvec(:,state+1)
        WRITE(iout,3) state+1, energy
     ELSE IF(i0stat == 'radial-packet') then
        CALL moment_data
        call rad_packet(wave_function)
     ELSE IF(i0stat == 'gaussian-pulse') THEN
        CALL moment_data
        CALL nr_packet(norm)
        CALL gauss_packet(wave_function,norm)
        CALL check_norm(wave_function,norm)
        CALL calc_moment(wave_function,t_calc)
     ELSE IF(i0stat == 'superpose') THEN
        CALL sppose
     ELSE
        CALL lnkerr('error in initial state')
     END IF
     call iosys('write real "initial state" to bec',nphy(1)*2,          &
                 wave_function,0,' ')
  ELSE
     WRITE(iout,2) t0
     CALL iosys ('read real solution from bec',nphy(1)*2,               &
                  wave_function,0,' ')
  END IF
  IF(prnton) THEN
     title='real and imaginary part of initial state'
     CALL print_psi(wave_function)
  END IF
  IF(t /= 1) THEN
     RETURN
  ELSE
     IF(imtime) THEN
        nfac=1.d0/ SQRT ( ddot(nphy(1),wave_function(1,1),1,            &
                                       wave_function(1,1),1)            &
                                        +                               &
                          ddot(nphy(1),wave_function(1,2),1,            &
                                       wave_function(1,2),1) )
        wave_function = nfac*wave_function
     END IF
  END IF
1    FORMAT(/,1X,'initial state construction at first time = ',e15.8,   &
    /,1X,'driver                                   = ',a24)
2    FORMAT(/,1X,'initial state from input file at time = ',e15.8)

3    FORMAT(/,5X,'initial state = ',i3,/,5X, 'energy        = ',e15.8)
END SUBROUTINE cp_psi_1d