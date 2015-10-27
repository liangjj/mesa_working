
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              MODULE initial_state
!**begin prologue     initial_state
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate the spatial part of the initial wavepacket.
!**references
!**                   
!**routines called    c_vect, rad_paket, gauss_paket, sppose
!**end prologue       initial_state
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              INTERFACE cp_psi
                    MODULE PROCEDURE cp_psi_1d_d, cp_psi_2d_d, cp_psi_3d_d
                              END INTERFACE cp_psi
!
                              INTERFACE c_vect                      
                    MODULE PROCEDURE c_vect_2d_d, c_vect_3d_d     
                              END INTERFACE c_vect
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck cp_psi_1d_d
!**begin prologue     cp_psi_1d_d
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
!**end prologue       cp_psi_1d_d
!
  SUBROUTINE cp_psi_1d_d(wave_function,t,t_calc)
  USE prop_global
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  USE normalize
  USE moment
  USE plot_wavefunction
  USE gaussian_wave_packet
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1))           :: wave_function
  INTEGER                              :: t, i 
  REAL*8                               :: norm, nfac
  REAL*8, DIMENSION(2)                 :: hfac
  LOGICAL                              :: dollar, logkey
  REAL*8                               :: ddot
  REAL*8                               :: t_calc
  CHARACTER (LEN=2)                    :: itoc
  CHARACTER (LEN=80)                   :: chrkey
  INTEGER                              :: intkey, state
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
     IF(i0stat.eq.'unperturbed-state-vector') then
        state=intkey(card,'initial-state',0,' ')
        energy=grid(1)%eigv_0(state+1)
        wave_function(:) = grid(1)%eigvec_0(:,state+1)   
        WRITE(iout,3) state+1, energy
     ELSE IF(i0stat == 'perturbed-state-vector') then
        state=intkey(card,'initial-state',0,' ')
        energy=grid(1)%eigv(state+1)
        wave_function(:) =   grid(1)%eigvec(:,state+1)
        CALL check_norm(wave_function,norm)
        WRITE(iout,3) state+1, energy
     ELSE IF(i0stat == 'radial-packet') then
        CALL moment_data
        call rad_packet(wave_function)
        CALL check_norm(wave_function,norm)
     ELSE IF(i0stat == 'gaussian-pulse') THEN
        CALL moment_data
        CALL nr_packet(norm)
        CALL gauss_packet(wave_function,norm)
        CALL check_norm(wave_function,norm)
     ELSE
        CALL lnkerr('error in initial state')
     END IF
     call iosys('write real "initial state" to bec',nphy(1),            &
                 wave_function,0,' ')
     call iosys('write real solution to bec',nphy(1),                   &
                 wave_function,0,' ')
  ELSE
     WRITE(iout,2) t0
     CALL iosys ('read real solution from bec',nphy(1),                 &
                  wave_function,0,' ')
  END IF
  IF(prnton) THEN
     title='initial state'
     CALL print_psi(wave_function)
  END IF
1 FORMAT(/,1X,'initial state construction at first time = ',e15.8,   &
         /,1X,'driver                                   = ',a24)
2 FORMAT(/,1X,'initial state from input file at time = ',e15.8)

3 FORMAT(/,5X,'initial state = ',i3,/,5X, 'energy        = ',e15.8)
END SUBROUTINE cp_psi_1d_d
!*deck cp_psi_2d_d
!**begin prologue     cp_psi_2d_d
!**date written       040706   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate the spatial part of the initial wavepacket.
!**references
!**                   
!**routines called    c_vect, rad_packet, gauss_packet, sppose
!**end prologue       cp_psi_2d_d
!
  SUBROUTINE cp_psi_2d_d(wave_function,t,t_calc)
  USE prop_global
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  USE normalize
  USE moment
  USE plot_wavefunction
  USE gaussian_wave_packet
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1))     :: wave_function
  INTEGER                                :: t, i 
  REAL*8                                 :: norm, nfac
  REAL*8, DIMENSION(2)                   :: hfac
  LOGICAL                                :: dollar, logkey
  REAL*8                                 :: ddot
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
     IF(i0stat.eq.'unperturbed-state-vector') then
        state=intkey(card,'initial-state',0,' ')
        ALLOCATE(itemp(nphy(2)*nphy(1),2),etemp(nphy(2)*nphy(1)))
        CALL c_vect(wave_function,etemp,itemp,state+1)
        CALL check_norm(wave_function,norm)
        DEALLOCATE(itemp,etemp)
     ELSE IF(i0stat == 'gaussian-pulse') THEN
        CALL moment_data
        CALL nr_packet(norm)
        CALL gauss_packet(wave_function,norm)
        CALL check_norm(wave_function,norm)
     ELSE
        CALL lnkerr('error in initial state')
     END IF
     call iosys('write real "initial state" to bec',nphy(2)*nphy(1),         &
                 wave_function,0,' ')
  ELSE
     WRITE(iout,2) t0
     CALL iosys ('read real solution from bec',nphy(2)*nphy(1),              &
                  wave_function,0,' ')
  END IF
  IF(prnton) THEN
     title='initial state'
     CALL print_psi(wave_function)
  END IF
1 FORMAT(/,1X,'initial state construction at first time = ',e15.8,        &
         /,1X,'driver                                   = ',a24)
2 FORMAT(/,1X,'initial state from input file at time = ',e15.8)
END SUBROUTINE cp_psi_2d_d
!*deck cp_psi_3d_d
!**begin prologue     cp_psi_3d_d
!**date written       040706   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate the spatial part of the initial wavepacket.
!**references
!**                   
!**routines called    c_vect, rad_paket, gauss_paket, sppose
!**end prologue       initial_state
!
  SUBROUTINE cp_psi_3d_d(wave_function,t,t_calc)
  USE prop_global
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  USE normalize
  USE moment
  USE plot_wavefunction
  USE gaussian_wave_packet
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))       :: wave_function
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))       :: scratch_vector
  INTEGER                                          :: t, i 
  REAL*8                                           :: norm, nfac
  REAL*8, DIMENSION(2)                             :: hfac
  LOGICAL                                          :: dollar, logkey
  REAL*8                                           :: ddot
  REAL*8                                           :: t_calc
  CHARACTER (LEN=2)                                :: itoc
  CHARACTER (LEN=80)                               :: chrkey
  INTEGER                                          :: intkey, state
  INTEGER                                          :: root
  REAL*8, DIMENSION(:),    ALLOCATABLE             :: etemp
  INTEGER, DIMENSION(:,:), ALLOCATABLE             :: itemp
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
     IF(i0stat.eq.'unperturbed-state-vector') then
        state=intkey(card,'initial-state',0,' ')
        ALLOCATE(itemp(nphy(3)*nphy(2)*nphy(1),3),                     &
                 etemp(nphy(3)*nphy(2)*nphy(1)))
        CALL c_vect(wave_function,etemp,itemp,state+1)
        CALL check_norm(wave_function,norm)
        DEALLOCATE(itemp,etemp)
     ELSE IF(i0stat == 'gaussian-pulse') THEN
        CALL moment_data
        CALL nr_packet(norm)
        CALL gauss_packet(wave_function,norm)
        CALL check_norm(wave_function,norm)
     ELSE
        CALL lnkerr('error in initial state')
     END IF
     call iosys('write real "initial state" to bec',                   &
                 nphy(3)*nphy(2)*nphy(1),wave_function,0,' ')
  ELSE
     WRITE(iout,2) t0
     CALL iosys ('read real solution from bec',                        &
                  nphy(3)*nphy(2)*nphy(1),wave_function,0,' ')
  END IF
  IF(prnton) THEN
     title='initial state'
     CALL print_psi(wave_function)
  END IF
1 FORMAT(/,1X,'initial state construction at first time = ',e15.8,  &
         /,1X,'driver                                   = ',a24)
2 FORMAT(/,1X,'initial state from input file at time = ',e15.8)

3 FORMAT(/,5X,'initial state = ',i3,/,5X, 'energy        = ',e15.8)
END SUBROUTINE cp_psi_3d_d
!deck c_vect_2d_d.f
!**begin prologue     c_vect_2d_d
!**date written       040607   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           hamiltonian
!**author             schneider, barry (nsf)
!**source             3-dim
!**purpose            Sort eigenvectors according energy and place
!***                  the desired root into initial state.
!**references
!**routines called
!**end prologue       c_vect_2d_d
  SUBROUTINE c_vect_2d_d(wave_function,etemp,itemp,root)
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                 :: root
  REAL*8, DIMENSION(nphy(2),nphy(1))      :: wave_function
  REAL*8, DIMENSION(nphy(2)*nphy(1))      :: etemp
  INTEGER, DIMENSION(nphy(2)*nphy(1),2)   :: itemp
  REAL*8                                  :: tmp
  INTEGER                                 :: i, j, k, i1, j1
  INTEGER                                 :: ii, count
!
!       The eig and ind arrays are destroyed by this routine.
!
!
!       Set up the index array
!
  count=0 
  do i=1,nphy(1)
     do j=1,nphy(2)
        count=count+1
        itemp(count,1)=i
        itemp(count,2)=j
     END DO
  END DO
!
!        Fill the eigenvalue array with all of the eigenvalues of the zeroth
!        order Hamiltonian.  These are simply the sum of the eigenvalues of
!        each dimension
!
  count=0
  DO  i=1,nphy(1)
      do j=1,nphy(2)
         count = count + 1
         etemp(count) = grid(1)%eigv_0(i) + grid(2)%eigv_0(j)
      END DO
  END DO
!
!        Sort the eigenvalues so that we can find the smallest
!
  DO  ii=2,n3d
      i=ii-1
      k=i
      tmp=etemp(i)
      i1=itemp(i,1)
      j1=itemp(i,2)
      DO  j=ii,n3d
          IF(etemp(j) < tmp) THEN
             k=j
             tmp=etemp(j)
          END IF
      END DO
      IF(k /= i) THEN
         itemp(i,1) = itemp(k,1)
         itemp(i,2) = itemp(k,2)
         itemp(k,1)=i1
         itemp(k,2)=j1
         etemp(k) = etemp(i)
         etemp(i) = tmp
      END IF
  END DO
!
!        Now that we know which eigenvector in each dimension fill the psi0
!        arrray with the appropriate eigenvector as the product.
!
  DO  i=1,nphy(1)
      DO j=1,nphy(2)
         wave_function(j,i) =  grid(1)%eigvec_0(i,itemp(root,1))      &
                                             *                        &
                               grid(2)%eigvec_0(j,itemp(root,2))
      END DO
  END DO
  energy=etemp(root)
  WRITE(iout,1) root, energy
1    FORMAT(/,5X,'initial state = ',i3,/,5X, 'energy        = ',e15.8)
END SUBROUTINE c_vect_2d_d
!deck c_vect_3d_d.f
!**begin prologue     c_vect_3d_d
!**date written       040607   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           hamiltonian
!**author             schneider, barry (nsf)
!**source             3-dim
!**purpose            Sort eigenvectors according energy and place
!***                  the desired root into initial state.
!**references
!**routines called
!**end prologue       c_vect_3d_d
  SUBROUTINE c_vect_3d_d(wave_function,etemp,itemp,root)
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                         :: root
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))      :: wave_function
  REAL*8, DIMENSION(nphy(3)*nphy(2)*nphy(1))      :: etemp
  INTEGER, DIMENSION(nphy(3)*nphy(2)*nphy(1),3)   :: itemp
  REAL*8                                          :: tmp
  INTEGER                                         :: i, j, k, i1, j1, k1
  INTEGER                                         :: ii, count
!
!       The eig and ind arrays are destroyed by this routine.
!
  count=0
  do i=1,nphy(1)
     do j=1,nphy(2)
        do k=1,nphy(3)
           count=count+1
           itemp(count,1)=i
           itemp(count,2)=j
           itemp(count,3)=k
        END DO
     END DO
  END DO 
  count=0
  DO  i=1,nphy(1)
      do j=1,nphy(2)
         do k=1,nphy(3)
            count = count + 1
            etemp(count) = grid(1)%eigv_0(i) +                         &
                           grid(2)%eigv_0(j) +                         &
                           grid(3)%eigv_0(k)
         END DO
      END DO
  END DO
  DO  ii=2,n3d
      i=ii-1
      k=i
      tmp=etemp(i)
      i1=itemp(i,1)
      j1=itemp(i,2)
      k1=itemp(i,3)
      DO  j=ii,n3d
          IF(etemp(j) < tmp) THEN
             k=j
             tmp=etemp(j)
          END IF
      END DO
      IF(k /= i) THEN
         itemp(i,1)=itemp(k,1)
         itemp(i,2)=itemp(k,2)
         itemp(i,3)=itemp(k,3)
         itemp(k,1)=i1
         itemp(k,2)=j1
         itemp(k,3)=k1
         etemp(k) = etemp(i)
         etemp(i) = tmp
     END IF
  END DO
  DO  i=1,nphy(1)
      do j=1,nphy(2)
         do k=1,nphy(3)
            wave_function(k,j,i) = grid(1)%eigvec_0(i,itemp(root,1)) * &
                                   grid(2)%eigvec_0(j,itemp(root,2)) * &
                                   grid(3)%eigvec_0(k,itemp(root,3))
         END DO
      END DO
  END DO
  energy=etemp(root)
  WRITE(iout,1) root, energy
1    FORMAT(/,5X,'initial state = ',i3,/,5X, 'energy        = ',e15.8)
END SUBROUTINE c_vect_3d_d
! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Subroutine for Gaussian Wavepacket Normalization}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck nr_packet.f
!**begin prologue     nr_packet
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate normailzation integrals for gaussian wavepacket
!**references
!**routines called
!**end prologue       nr_packet
  SUBROUTINE nr_packet(norm)
  USE arnoldi_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8                                 :: norm
  REAL*8, DIMENSION(3)                   :: ov 
  INTEGER                                :: i
  IF(system == 'cartesian') THEN
     norm=1.d0
     DO i=1,spdim
        call ov1_quad(ov(i),grid(i)%pt,grid(i)%wt, &
                      x_0(i),alpha(i),sigma(i),powr(i),nphy(i))
        norm=norm*ov(i)
     END DO
  ELSE
     call ov2_quad(norm,grid(1)%pt,grid(1)%wt, &
                   x_0(1),alpha(1),sigma(1),nphy(1))    
  END IF
  norm=1.d0/SQRT(norm)
  WRITE(iout,1) norm
1    FORMAT(/,1X,'overlap integral for gaussian packet = ',e15.8)
END SUBROUTINE nr_packet
END MODULE initial_state
