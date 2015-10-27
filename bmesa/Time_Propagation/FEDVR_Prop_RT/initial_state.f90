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
                      INTERFACE cp_psi
           MODULE PROCEDURE cp_psi_1d_z, cp_psi_2d_z, cp_psi_3d_z
                  END INTERFACE cp_psi
                      INTERFACE c_vect
           MODULE PROCEDURE c_vect_2d_z, c_vect_3d_z
                  END INTERFACE c_vect
                      INTERFACE sppose
           MODULE PROCEDURE sppose_1d_z, sppose_2d_z, sppose_3d_z
                  END INTERFACE sppose
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck cp_psi_1d_z
!**begin prologue     cp_psi_1d_z
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
!**end prologue       cp_psi_1d_z
!
  SUBROUTINE cp_psi_1d_z(wave_function,scratch_vector,t,t_calc)
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
        energy=grid(1)%eigv_0(state+1)
        wave_function(:,1) = grid(1)%eigvec_0(:,state+1)   
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
        wave_function(:,1) = grid(1)%eigvec(:,state+1)   
        wave_function(:,2) = - sin(energy*t0*hfac(1))                   &
                                         *                              &
                               wave_function(:,1)
        wave_function(:,1) =   cos(energy*t0*hfac(1))                   &
                                         *                              &
                               wave_function(:,1)
        CALL check_norm(wave_function,norm)
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
        CALL sppose(wave_function)
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
END SUBROUTINE cp_psi_1d_z
!*deck cp_psi_2d_z
!**begin prologue     cp_psi_2d_z
!**date written       040706   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate the spatial part of the initial wavepacket.
!**references
!**                   
!**routines called    c_vect, rad_packet, gauss_packet, sppose
!**end prologue       cp_psi_2d_z
!
  SUBROUTINE cp_psi_2d_z(wave_function,scratch_vector,t,t_calc)
  USE prop_global
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  USE normalize
  USE moment
  USE plot_wavefunction
  USE gaussian_wave_packet
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1),2)   :: wave_function
  REAL*8, DIMENSION(nphy(2),nphy(1),2)   :: scratch_vector
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
     IF(i0stat == 'one') THEN
        wave_function(:,:,1) =    cos (t0*t0*hfac(2))
        wave_function(:,:,2) =  - sin (t0*t0*hfac(2))
!        If the right hand side is a state vector, we need some scratch
!        memory to sort the eigenstates by energy.  Go get it.
     ELSE IF(i0stat.eq.'unperturbed-state-vector') then
        state=intkey(card,'initial-state',0,' ')
        ALLOCATE(itemp(nphy(2)*nphy(1),2),etemp(nphy(2)*nphy(1)))
        CALL c_vect(wave_function(:,:,1),etemp,itemp,state+1)
!   
!       Use scratch_vector as a temporary
!          
        scratch_vector(:,:,1) =   wave_function(:,:,1)                  &
                                       *                                &
                                  cos(energy*t0*hfac(1))   
        scratch_vector(:,:,2) = - wave_function(:,:,1)                  &
                                       *                                &
                                  sin(energy*t0*hfac(1)) 
!
!       Now store it
! 
        wave_function = scratch_vector
        DEALLOCATE(itemp,etemp)
     ELSE IF(i0stat == 'radial-packet') then
        CALL moment_data
        call rad_paket(wave_function)
     ELSE IF(i0stat == 'gaussian-pulse') THEN
        CALL moment_data
        CALL nr_packet(norm)
        CALL gauss_packet(wave_function,norm)
        CALL check_norm(wave_function,norm,f_1)
        CALL calc_moment(wave_function,t_calc,f_2)
     ELSE IF(i0stat == 'superpose') THEN
        CALL sppose(wave_function)
     ELSE
        CALL lnkerr('error in initial state')
     END IF
     call iosys('write real "initial state" to bec',nphy(2)*nphy(1)*2,       &
                 wave_function,0,' ')
  ELSE
     WRITE(iout,2) t0
     CALL iosys ('read real solution from bec',nphy(2)*nphy(1)*2,            &
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
        nfac=1.d0/ SQRT ( ddot(nphy(2)*nphy(1),wave_function(1,1,1),         &
                                               1,                            &
                                               wave_function(1,1,1),         &
                                               1)                            &
                                               +                             &
                          ddot(nphy(2)*nphy(1),wave_function(1,1,2),         &
                                               1,                            &
                                               wave_function(1,1,2),         &
                                               1) )
        wave_function = nfac*wave_function
     END IF
  END IF
1    FORMAT(/,1X,'initial state construction at first time = ',e15.8,        &
    /,1X,'driver                                   = ',a24)
2    FORMAT(/,1X,'initial state from input file at time = ',e15.8)
END SUBROUTINE cp_psi_2d_z
!*deck cp_psi_3d_z
!**begin prologue     cp_psi_3d_z
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
  SUBROUTINE cp_psi_3d_z(wave_function,scratch_vector,t,t_calc)
  USE prop_global
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  USE normalize
  USE moment
  USE plot_wavefunction
  USE gaussian_wave_packet
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)     :: wave_function
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)     :: scratch_vector
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
     IF(i0stat == 'one') THEN
        wave_function(:,:,:,1) =    cos (t0*t0*hfac(2))
        wave_function(:,:,:,2) =  - sin (t0*t0*hfac(2))
!        soln = wave_function
!        If the right hand side is a state vector, we need some scratch
!        memory to sort the eigenstates by energy.  Go get it.
     ELSE IF(i0stat.eq.'unperturbed-state-vector') then
        state=intkey(card,'initial-state',0,' ')
        ALLOCATE(itemp(nphy(3)*nphy(2)*nphy(1),3),                     &
                 etemp(nphy(3)*nphy(2)*nphy(1)))
         CALL c_vect(wave_function(:,:,:,1),etemp,itemp,state+1)
!   
!       Use scratch_vector as a temporary
!          
        scratch_vector(:,:,:,1) =   wave_function(:,:,:,1)             &
                                         *                             &
                                    cos(energy*t0*hfac(1))   
        scratch_vector(:,:,:,2) = - wave_function(:,:,:,1)             &
                                         *                             &
                                    sin(energy*t0*hfac(1)) 
!
!       Now store it
! 
        wave_function = scratch_vector
        DEALLOCATE(itemp,etemp)
     ELSE IF(i0stat == 'radial-packet') then
        CALL moment_data
        call rad_packet(wave_function)
     ELSE IF(i0stat == 'gaussian-pulse') THEN
        CALL moment_data
        CALL nr_packet(norm)
        CALL gauss_packet(wave_function,norm)
        CALL check_norm(wave_function,norm,f_1,f_2)
        CALL calc_moment(wave_function,t_calc,f_2,f_3)
     ELSE IF(i0stat == 'superpose') THEN
        CALL sppose(wave_function)
     ELSE
        CALL lnkerr('error in initial state')
     END IF
     call iosys('write real "initial state" to bec',                   &
                 nphy(3)*nphy(2)*nphy(1)*2,wave_function,0,' ')
  ELSE
     WRITE(iout,2) t0
     CALL iosys ('read real solution from bec',                        &
                  nphy(3)*nphy(2)*nphy(1)*2,wave_function,0,' ')
  END IF
  IF(prnton) THEN
     title='real and imaginary part of initial state'
     CALL print_psi(wave_function)
  END IF
  IF(t /= 1) THEN
     RETURN
  ELSE
     IF(imtime) THEN
        nfac=1.d0/ SQRT ( ddot(nphy(3)*nphy(2)*nphy(1),                &
                                       wave_function(1,1,1,1),         &
                                       1,                              &
                                       wave_function(1,1,1,1),         &
                                       1)                              &
                                              +                        &
                          ddot(nphy(3)*nphy(2)*nphy(1),                &
                                       wave_function(1,1,1,2),         &
                                       1,                              &
                                       wave_function(1,1,1,2),         &
                                       1) )
        wave_function = nfac*wave_function
     END IF
  END IF
1    FORMAT(/,1X,'initial state construction at first time = ',e15.8,  &
    /,1X,'driver                                   = ',a24)
2    FORMAT(/,1X,'initial state from input file at time = ',e15.8)

3    FORMAT(/,5X,'initial state = ',i3,/,5X, 'energy        = ',e15.8)
END SUBROUTINE cp_psi_3d_z
!deck c_vect_2d_z.f
!**begin prologue     c_vect_2d_z
!**date written       040607   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           hamiltonian
!**author             schneider, barry (nsf)
!**source             3-dim
!**purpose            Sort eigenvectors according energy and place
!***                  the desired root into initial state.
!**references
!**routines called
!**end prologue       c_vect_2d_z
  SUBROUTINE c_vect_2d_z(wave_function,etemp,itemp,root)
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
END SUBROUTINE c_vect_2d_z
!deck c_vect_3d_z.f
!**begin prologue     c_vect_3d_z
!**date written       040607   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           hamiltonian
!**author             schneider, barry (nsf)
!**source             3-dim
!**purpose            Sort eigenvectors according energy and place
!***                  the desired root into initial state.
!**references
!**routines called
!**end prologue       c_vect_3d_z_z
  SUBROUTINE c_vect_3d_z(wave_function,etemp,itemp,root)
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
END SUBROUTINE c_vect_3d_z
!deck sppose_1d_z.f
!***begin prologue     sppose
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            calculate zero time gaussian wavepacket
!***
!***references
!***routines called
!***end prologue       sppose_1d_z
  SUBROUTINE sppose_1d_z(wave_function)
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  LOGICAL                                    :: dollar
  INTEGER                                    :: ntot
  INTEGER, DIMENSION(3)                      :: n_val
  INTEGER                                    :: i, intkey
  REAL*8, DIMENSION(nphy(1),2)               :: wave_function
  ALLOCATE(c_loc(1))
  ALLOCATE(c_loc(1)%c(nphy(1)),c_loc(1)%phi(nphy(1)),c_loc(1)%l(nphy(1)))
  wave_function = 0.d0
  IF ( dollar('$states',card,title,inp) ) THEN
       IF(system == 'cartesian') THEN
           ntot=1
           n_val(1)=intkey(card,'number-of-x-superposed-states',1,' ')
           ntot=ntot*n_val(1)
           CALL intarr(card,'x-state-list',c_loc(1)%l,n_val(1),' ')
           CALL fparr(card,'x-state-coefficients',c_loc(1)%c,n_val(1),' ')
           CALL mk_phi(c_loc(1)%phi,grid(1)%eigvec_0,c_loc(1)%c,c_loc(1)%l, &
                       n_val(1),nphy(1),'x')
       ELSE
           n_val(1)=intkey(card,'number-of-r-superposed-states',1,' ')
           ntot=n_val(1)
           CALL intarr(card,'r-state-list',c_loc(1)%l,n_val(1),' ')
           CALL fparr(card,'r-state-coefficients',c_loc(1)%c,n_val(1),' ')
           CALL mk_phi(c_loc(1)%phi,grid(1)%eigvec_0,c_loc(1)%c,c_loc(1)%l, &
                       n_val(1),nphy(1),'r')
       END IF
  END IF
  wave_function(:,1) = wave_function(:,1) + c_loc(1)%phi(:)
  DEALLOCATE(c_loc(1)%c,c_loc(1)%phi,c_loc(1)%l)
  DEALLOCATE(c_loc)
END SUBROUTINE sppose_1d_z
!deck sppose_2d_z.f
!***begin prologue     sppose_2d_z
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            calculate zero time gaussian wavepacket
!***
!***references
!***routines called
!***end prologue       sppose_2d_z
  SUBROUTINE sppose_2d_z(wave_function)
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  LOGICAL                                        :: dollar
  INTEGER                                        :: ntot
  INTEGER, DIMENSION(3)                          :: n_val
  INTEGER                                        :: i, j, intkey
  REAL*8, DIMENSION(nphy(2),nphy(1),2)           :: wave_function
  ALLOCATE(c_loc(2))
  DO i=1,2
     ALLOCATE(c_loc(i)%c(nphy(i)),c_loc(i)%phi(nphy(i)),c_loc(i)%l(nphy(i)))
  END DO
  wave_function = 0.d0
  IF ( dollar('$states',card,title,inp) ) THEN
       ntot=1
       n_val(1)=intkey(card,'number-of-x-superposed-states',1,' ')
       ntot=ntot*n_val(1)
       CALL intarr(card,'x-state-list',c_loc(1)%l,n_val(1),' ')
       CALL fparr(card,'x-state-coefficients',c_loc(1)%c,n_val(1),' ')
       CALL mk_phi(c_loc(1)%phi,grid(1)%eigvec_0,c_loc(1)%c,c_loc(1)%l, &
                   n_val(1),nphy(1),'x')
       n_val(2)=intkey(card,'number-of-y-superposed-states',1,' ')
       ntot=ntot*n_val(2)
       CALL intarr(card,'y-state-list',c_loc(2)%l,n_val(2),' ')
       CALL fparr(card,'y-state-coefficients',c_loc(2)%c,n_val(2),' ')
       CALL mk_phi(c_loc(2)%phi,grid(2)%eigvec_0,c_loc(2)%c,c_loc(2)%l, &
                   n_val(2),nphy(2),'y')
       DO i=1,nphy(1)
          DO j=1,nphy(2)
             wave_function(j,i,1) = wave_function(j,i,1) +              &
                                    c_loc(2)%phi(j) +c_loc(1)%phi(i)
          END DO
       END DO
  END IF
  DO i=1,2
     DEALLOCATE(c_loc(i)%c,c_loc(i)%phi,c_loc(i)%l)
  END DO
  DEALLOCATE(c_loc)
END SUBROUTINE sppose_2d_z
!deck sppose_3d.f
!***begin prologue     sppose_3d_z
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            calculate zero time gaussian wavepacket
!***
!***references
!***routines called
!***end prologue       sppose_3d_z
  SUBROUTINE sppose_3d_z(wave_function)
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  LOGICAL                                        :: dollar
  INTEGER                                        :: ntot
  INTEGER, DIMENSION(3)                          :: n_val
  INTEGER                                        :: i, j, k, intkey
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: wave_function
  ALLOCATE(c_loc(2))
  DO i=1,3
     ALLOCATE(c_loc(i)%c(nphy(i)),c_loc(i)%phi(nphy(i)),c_loc(i)%l(nphy(i)))
  END DO
  wave_function = 0.d0
  IF ( dollar('$states',card,title,inp) ) THEN
       ntot=1
       n_val(1)=intkey(card,'number-of-x-superposed-states',1,' ')
       ntot=ntot*n_val(1)
       CALL intarr(card,'x-state-list',c_loc(1)%l,n_val(1),' ')
       CALL fparr(card,'x-state-coefficients',c_loc(1)%c,n_val(1),' ')
       CALL mk_phi(c_loc(1)%phi,grid(1)%eigvec_0,c_loc(1)%c,c_loc(1)%l, &
                   n_val(1),nphy(1),'x')
       n_val(2)=intkey(card,'number-of-y-superposed-states',1,' ')
       ntot=ntot*n_val(2)
       CALL intarr(card,'y-state-list',c_loc(2)%l,n_val(2),' ')
       CALL fparr(card,'y-state-coefficients',c_loc(2)%c,n_val(2),' ')
       CALL mk_phi(c_loc(2)%phi,grid(2)%eigvec_0,c_loc(2)%c,c_loc(2)%l, &
                   n_val(2),nphy(2),'y')
       n_val(3)=intkey(card,'number-of-z-superposed-states',1,' ')
       ntot=ntot*n_val(3)
       CALL intarr(card,'z-state-list',c_loc(3)%l,n_val(3),' ')
       CALL fparr(card,'z-state-coefficients',c_loc(3)%c,n_val(3),' ')
       CALL mk_phi(c_loc(3)%phi,grid(3)%eigvec_0,c_loc(3)%c,c_loc(3)%l, &
                   n_val(3),nphy(3),'z')
       DO i=1,nphy(1)
          DO j=1,nphy(2)
             Do k=1,nphy(3)
                wave_function(k,j,i,1) = wave_function(k,j,i,1)            &
                                               +                           &
                                         c_loc(3)%phi(k)                   &
                                               +                           &
                                         c_loc(2)%phi(j)                   &
                                               +                           &
                                         c_loc(1)%phi(i)
             END DO 
          END DO
       END DO
  END IF
  DO i=1,3
     DEALLOCATE(c_loc(i)%c,c_loc(i)%phi,c_loc(i)%l)
  END DO
  DEALLOCATE(c_loc)
END SUBROUTINE sppose_3d_z
END MODULE initial_state
