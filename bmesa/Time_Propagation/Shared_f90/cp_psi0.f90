! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Cp_psi0}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck cp_psi0.f
!**begin prologue     cp_psi0
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate the spatial part of the initial wavepacket.
!**references
!**                   
!**routines called
!**end prologue       cp_psi0
  SUBROUTINE cp_psi0(t)
  USE arnoldi_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: t 
  REAL*8                                 :: f1, f2, fpkey
  REAL*8                                 :: norm, st, nfac
  LOGICAL                                :: dollar, logkey
  COMPLEX*16                             :: cdotc
  CHARACTER (LEN=2)                      :: itoc
  CHARACTER (LEN=80)                     :: chrkey
  INTEGER                                :: intkey, state
  REAL*8, DIMENSION(:),    ALLOCATABLE   :: etemp
  INTEGER, DIMENSION(:,:), ALLOCATABLE   :: itemp
  SAVE norm
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
        psi0 = EXP (-eye*t0*t0*.5D0)
        soln_0 = psi0
!        If the right hand side is a state vector, we need some scratch
!        memory to sort the eigenstates by energy.  Go get it.
     ELSE IF(i0stat.eq.'unperturbed-state-vector') then
        state=intkey(card,'initial-state',0,' ')
        IF(spdim > 1) THEN
           ALLOCATE(itemp(n3d,spdim),etemp(n3d))
        END IF
        CALL c_vect0(etemp,itemp,state+1)
        psi0 = psi0*EXP(-eye*energy*t0)
        soln_0 = psi0
        IF(spdim > 1) THEN
           DEALLOCATE(itemp,etemp)
        END IF
     ELSE IF(i0stat == 'perturbed-state-vector') then
        state=intkey(card,'initial-state',0,' ')
        energy=grid(1)%eigv(state+1)
        psi0 = EXP(-eye*energy*t0)*grid(1)%eigvec(:,state+1)
        soln_0 = psi0    
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
        soln_0 = psi0
     ELSE IF(i0stat == 'gaussian-pulse') THEN
        CALL fparr(card,'alpha',alpha,3,' ')
        CALL fparr(card,'sigma',sigma,3,' ')
        CALL fparr(card,'x_0',x_0,3,' ')
        CALL fparr(card,'beta',beta,3,' ')
        CALL nr_paket(norm)
        CALL gauss_paket(norm)
        soln_0 = psi0
     ELSE IF(i0stat == 'superpose') THEN
        CALL sppose
        soln_0 = psi0
     ELSE
        CALL lnkerr('error in initial state')
     END IF
  ELSE
     WRITE(iout,2) t0
     CALL iosys ('read real solution from bec',n3d*2,psi0,0,' ')
  END IF
  IF(prnton) THEN
     title='initial state'
     CALL prntcmn(title,psi0,n3d,1,n3d,1,iout,'e')
  END IF
  IF(t /= 1) THEN
     RETURN
  ELSE
     IF(imtime) THEN
        nfac=1.d0/ SQRT ( cdotc(n3d,psi0,1,psi0,1) )
        psi0 = nfac*psi0
     END IF
  END IF
1    FORMAT(/,1X,'initial state construction at first time = ',e15.8,  &
    /,1X,'driver                                   = ',a24)
2    FORMAT(/,1X,'initial state from input file at time = ',e15.8)

3    FORMAT(/,5X,'initial state = ',i3,/,5X, 'energy        = ',e15.8)
END SUBROUTINE cp_psi0
