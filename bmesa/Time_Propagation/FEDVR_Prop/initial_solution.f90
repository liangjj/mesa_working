!deck initial_solution.f
!**begin prologue     initial_solution
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate the spatial part of the initial wavepacket.
!**references
!**                   
!**routines called    c_vect_08, rad_paket_08, gauss_paket_08, sppose_08
!**end prologue       initial_solution
  SUBROUTINE initial_solution(t,psi,v_scr)
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: t 
  REAL*8, DIMENSION(n3d,2)               :: psi, v_scr 
  INTEGER                                :: i 
  REAL*8                                 :: f1, f2, fpkey
  REAL*8                                 :: norm, st, nfac
  REAL*8, DIMENSION(2)                   :: hfac
  LOGICAL                                :: dollar, logkey
  REAL*8                                 :: sdot
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
!        soln_08 = psi_08
!        If the right hand side is a state vector, we need some scratch
!        memory to sort the eigenstates by energy.  Go get it.
     ELSE IF(i0stat.eq.'unperturbed-state-vector') then
        state=intkey(card,'initial-state',0,' ')
        IF(spdim > 1) THEN
           ALLOCATE(itemp(n3d,spdim),etemp(n3d))
        END IF
        CALL c_vect(psi,etemp,itemp,state+1)
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
        WRITE(iout,3) state+1, energy
     ELSE IF(i0stat == 'radial-packet') then
        alpha(1)=fpkey(card,'alpha',1.d0,' ')
        sigma(1)=fpkey(card,'sigma',1.d0,' ')
        x_0(1)=fpkey(card,'x_0',0.d0,' ')
        beta(1)=fpkey(card,'beta',0.d0,' ')
        powr=intkey(card,'power',1,' ')
        typ_pak=chrkey(card,'type-radial-packet', &
                            'exponential',' ')
        call rad_paket
     ELSE IF(i0stat == 'gaussian-pulse') THEN
        CALL fparr(card,'alpha',alpha,3,' ')
        CALL fparr(card,'sigma',sigma,3,' ')
        CALL fparr(card,'x_0',x_0,3,' ')
        CALL fparr(card,'beta',beta,3,' ')
        CALL nr_paket(norm)
        CALL gauss_paket_08(norm)
        IF(spdim==1) THEN
           CALL chk_nrm_08_1d(psi_08)
        ELSE IF(spdim==2) THEN
           CALL chk_nrm_08_2d(psi_08,fac)
        ELSE IF(spdim==3) THEN
           CALL chk_nrm_08_3d(psi_08,fac,v_1)
        END IF
        CALL moment_08
!        soln_08 = psi_08
     ELSE IF(i0stat == 'superpose') THEN
        CALL sppose_08
!        soln_08 = psi_08
     ELSE
        CALL lnkerr('error in initial state')
     END IF
     call iosys('write real "initial state" to bec',n3d*2,psi_08,0,' ')
  ELSE
     WRITE(iout,2) t0
     CALL iosys ('read real solution from bec',n3d*2,psi_08,0,' ')
  END IF
  IF(prnton) THEN
     title='real and imaginary part of initial state'
     CALL plot_wavefunction(psi_08)
  END IF
  IF(t /= 1) THEN
     RETURN
  ELSE
     IF(imtime) THEN
        nfac=1.d0/ SQRT ( sdot(n3d,psi_08(1,1),1,psi_08(1,1),1) + &
                          sdot(n3d,psi_08(1,2),1,psi_08(1,2),1) )
        psi_08 = nfac*psi_08
     END IF
  END IF
1    FORMAT(/,1X,'initial state construction at first time = ',e15.8,  &
    /,1X,'driver                                   = ',a24)
2    FORMAT(/,1X,'initial state from input file at time = ',e15.8)

3    FORMAT(/,5X,'initial state = ',i3,/,5X, 'energy        = ',e15.8)
END SUBROUTINE initial_solution
