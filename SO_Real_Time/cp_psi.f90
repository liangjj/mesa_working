!deck cp_psi.f
!**begin prologue     cp_psi
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate the spatial part of the initial wavepacket.
!**references
!**                   
!**routines called    c_vect, rad_paket, gauss_paket, sppose
!**end prologue       cp_psi
  SUBROUTINE cp_psi(t,t_calc)
  USE prop_global
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  USE normalize
  USE moment
  USE plot_wavefunction
  IMPLICIT NONE
  INTEGER                                :: t, i 
  REAL*8                                 :: f1, f2, fpkey
  REAL*8                                 :: norm, st, nfac
  REAL*8, DIMENSION(2)                   :: hfac
  LOGICAL                                :: dollar, logkey
  REAL*8                                 :: sdot
  REAL*8                                 :: t_calc
  CHARACTER (LEN=2)                      :: itoc
  CHARACTER (LEN=80)                     :: chrkey
  INTEGER                                :: intkey, state
  REAL*8, DIMENSION(:),    ALLOCATABLE   :: etemp
  INTEGER, DIMENSION(:,:), ALLOCATABLE   :: itemp
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
        psi(:,1) =    cos (t0*t0*hfac(2))
        psi(:,2) =  - sin (t0*t0*hfac(2))
!        soln = psi
!        If the right hand side is a state vector, we need some scratch
!        memory to sort the eigenstates by energy.  Go get it.
     ELSE IF(i0stat.eq.'unperturbed-state-vector') then
        state=intkey(card,'initial-state',0,' ')
        IF(spdim > 1) THEN
           ALLOCATE(itemp(n3d,spdim),etemp(n3d))
        END IF
        CALL c_vect(etemp,itemp,state+1)
!   
!       Use v_scr as a temporary
!          
        v_scr(:,1) =   psi(:,1) * cos(energy*t0*hfac(1))   
        v_scr(:,2) = - psi(:,1) * sin(energy*t0*hfac(1)) 
!
!       Now store it
! 
        psi = v_scr
        IF(spdim > 1) THEN
           DEALLOCATE(itemp,etemp)
        END IF
     ELSE IF(i0stat == 'perturbed-state-vector') then
        state=intkey(card,'initial-state',0,' ')
        energy=grid(1)%eigv(state+1)
        psi(:,1) =   cos(energy*t0*hfac(1))*grid(1)%eigvec(:,state+1)
        psi(:,2) = - sin(energy*t0*hfac(1))*grid(1)%eigvec(:,state+1)
!        soln = psi
        WRITE(iout,3) state+1, energy
     ELSE IF(i0stat == 'radial-packet') then
        CALL moment_data
        call rad_paket
     ELSE IF(i0stat == 'gaussian-pulse') THEN
        CALL moment_data
        CALL nr_paket(norm)
        CALL gauss_paket(norm)
        IF(spdim==1) THEN
           CALL check_norm_1d(psi)
           CALL moment_1d(psi,t_calc)
        ELSE IF(spdim==2) THEN
           CALL check_norm_2d(psi,fac)
           CALL moment_2d(psi,fac,t_calc)
        ELSE IF(spdim==3) THEN
           CALL check_norm_3d(psi,fac,v_1)
           CALL moment_3d(psi,fac,v_1,t_calc)
        END IF

!        soln = psi
     ELSE IF(i0stat == 'superpose') THEN
        CALL sppose
!        soln = psi
     ELSE
        CALL lnkerr('error in initial state')
     END IF
     call iosys('write real "initial state" to bec',n3d*2,psi,0,' ')
  ELSE
     WRITE(iout,2) t0
     CALL iosys ('read real solution from bec',n3d*2,psi,0,' ')
  END IF
  IF(prnton) THEN
     title='real and imaginary part of initial state'
     CALL print_wavefunction(psi)
  END IF
  IF(t /= 1) THEN
     RETURN
  ELSE
     IF(imtime) THEN
        nfac=1.d0/ SQRT ( sdot(n3d,psi(1,1),1,psi(1,1),1) + &
                          sdot(n3d,psi(1,2),1,psi(1,2),1) )
        psi = nfac*psi
     END IF
  END IF
1    FORMAT(/,1X,'initial state construction at first time = ',e15.8,  &
    /,1X,'driver                                   = ',a24)
2    FORMAT(/,1X,'initial state from input file at time = ',e15.8)

3    FORMAT(/,5X,'initial state = ',i3,/,5X, 'energy        = ',e15.8)
END SUBROUTINE cp_psi
