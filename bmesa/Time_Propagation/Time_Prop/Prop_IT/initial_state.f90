!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              MODULE initial_state
                              USE prop_global
                              USE dvrprop_global
                              USE dvr_shared
                              USE dvr_global
                              USE normalize
                              USE moment
                              USE plot_wavefunction
                              USE gaussian_wave_packet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
                              INTERFACE initial_vector
                    MODULE PROCEDURE initial_vector_d,                  &
                                     initial_vector_z
                              END INTERFACE initial_vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              INTERFACE sppose                      
                    MODULE PROCEDURE sppose_d,                          &
                                     sppose_z     
                              END INTERFACE sppose
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck initial_vector_d
!**begin prologue     initial_vector_d
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
  SUBROUTINE initial_vector_d(wave_function,t,t_calc)
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d)                           :: wave_function
  REAL*8, DIMENSION(n3d)                           :: scratch_vector
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
     IF(i0stat.eq.'unperturbed-state-vector') THEN
        state=intkey(card,'initial-state',0,' ')
        IF( spdim == 1) THEN
            energy=grid(1)%eigv_0(state+1)
            wave_function(:) = grid(1)%eigvec_0(:,state+1)   
        ELSE IF (spdim == 2 ) THEN
            ALLOCATE(itemp(n3d,2),etemp(n3d))
            CALL c_vect_2d_d(wave_function,etemp,itemp,state+1)
            DEALLOCATE(itemp,etemp)
        ELSE IF (spdim == 3) THEN
            ALLOCATE(itemp(n3d,3),etemp(n3d))
            CALL c_vect_3d_d(wave_function,etemp,itemp,state+1)
            DEALLOCATE(itemp,etemp)
        END IF
        CALL check_norm(wave_function,norm)
     ELSE IF(i0stat.eq.'perturbed-state-vector') THEN
        state=intkey(card,'initial-state',0,' ')
        IF( spdim == 1) THEN
            energy=grid(1)%eigv(state+1)
            wave_function(:) = grid(1)%eigvec(:,state+1)   
        END IF 
     ELSE IF(i0stat == 'gaussian-pulse') THEN
        CALL moment_data
        CALL nr_packet(norm)
        CALL gauss_packet(wave_function,norm)
        CALL check_norm(wave_function,norm)
     ELSE IF(i0stat == 'superpose' ) THEN
        CALL sppose(wave_function)
     ELSE
        CALL lnkerr('error in initial state')
     END IF
     call iosys('write real "initial state" to bec',n3d,wave_function,0,' ')
     call iosys('write real solution to bec',n3d,wave_function,0,' ')
  ELSE
     CALL iosys ('read real solution from bec',n3d,wave_function,0,' ')
  END IF
  IF(prnton) THEN
     title='initial state'
     CALL print_psi(wave_function)
  END IF
1 FORMAT(/,1X,'initial state construction at first time = ',e15.8,  &
         /,1X,'driver                                   = ',a24)
END SUBROUTINE initial_vector_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck initial_vector_z
!**begin prologue     initial_vector_z
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
  SUBROUTINE initial_vector_z(wave_function,t,t_calc)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n3d)                       :: wave_function
  COMPLEX*16, DIMENSION(n3d)                       :: scratch_vector
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
        IF (spdim == 2 ) THEN
            ALLOCATE(itemp(n3d,2),etemp(n3d))
            CALL c_vect_2d_z(wave_function,etemp,itemp,state+1)
            DEALLOCATE(itemp,etemp)
        ELSE IF (spdim == 3) THEN
            ALLOCATE(itemp(n3d,3),etemp(n3d))
            CALL c_vect_3d_z(wave_function,etemp,itemp,state+1)
            DEALLOCATE(itemp,etemp)
        END IF
        CALL check_norm(wave_function,norm)
     ELSE IF(i0stat == 'gaussian-pulse') THEN
        CALL moment_data
        CALL nr_packet(norm)
        CALL gauss_packet(wave_function,norm)
        CALL check_norm(wave_function,norm)
     ELSE IF(i0stat == 'superpose' ) THEN
        CALL sppose(wave_function)
     ELSE
        CALL lnkerr('error in initial state')
     END IF
     call iosys('write real "initial state" to bec',2*n3d,          &
                 wave_function,0,' ')
  ELSE
     CALL iosys ('read real solution from bec',2*n3d,               &
                  wave_function,0,' ')
  END IF
  IF(prnton) THEN
     title='initial state'
     CALL print_psi(wave_function)
  END IF
1 FORMAT(/,1X,'initial state construction at first time = ',e15.8,  &
         /,1X,'driver                                   = ',a24)
END SUBROUTINE initial_vector_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  IMPLICIT NONE
  INTEGER                                 :: root
  COMPLEX*16, DIMENSION(nphy(2),nphy(1))  :: wave_function
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!**end prologue       c_vect_3d_z
  SUBROUTINE c_vect_3d_z(wave_function,etemp,itemp,root)
  IMPLICIT NONE
  INTEGER                                         :: root
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1))  :: wave_function
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck sppose_d.f
!***begin prologue     sppose_d
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            calculate zero time gaussian wavepacket
!***
!***references
!***routines called
!***end prologue       sppose_d
  SUBROUTINE sppose_d(wave_function)
  IMPLICIT NONE
  LOGICAL                                        :: dollar
  INTEGER                                        :: ntot
  INTEGER, DIMENSION(3)                          :: n_val
  INTEGER                                        :: i, j, k, intkey
  REAL*8, DIMENSION(n3d)                         :: wave_function
  ALLOCATE(c_loc(spdim))
  DO i=1,spdim
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
       IF (spdim >= 2) THEN
           n_val(2)=intkey(card,'number-of-y-superposed-states',1,' ')
           ntot=ntot*n_val(2)
           CALL intarr(card,'y-state-list',c_loc(2)%l,n_val(2),' ')
           CALL fparr(card,'y-state-coefficients',c_loc(2)%c,n_val(2),' ')
           CALL mk_phi(c_loc(2)%phi,grid(2)%eigvec_0,c_loc(2)%c,c_loc(2)%l, &
                       n_val(2),nphy(2),'y')
       END IF
       IF (spdim >= 3) THEN
           n_val(3)=intkey(card,'number-of-z-superposed-states',1,' ')
           ntot=ntot*n_val(3)
           CALL intarr(card,'z-state-list',c_loc(3)%l,n_val(3),' ')
           CALL fparr(card,'z-state-coefficients',c_loc(3)%c,n_val(3),' ')
           CALL mk_phi(c_loc(3)%phi,grid(3)%eigvec_0,c_loc(3)%c,c_loc(3)%l, &
                       n_val(3),nphy(3),'z')
       END IF
  END IF
  IF (spdim ==1 ) THEN
      CALL wfn_fil_1d_dd(wave_function,c_loc(1)%phi)
  ELSE IF( spdim == 2) THEN
      CALL wfn_fil_2d_dd(wave_function,c_loc(1)%phi,c_loc(2)%phi)
  ELSE IF( spdim == 3) THEN
      CALL wfn_fil_3d_dd(wave_function,c_loc(1)%phi,c_loc(2)%phi,c_loc(3)%phi)
  END IF
  DO i=1,spdim
     DEALLOCATE(c_loc(i)%c,c_loc(i)%phi,c_loc(i)%l)
  END DO
  DEALLOCATE(c_loc)
END SUBROUTINE sppose_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck sppose_z.f
!***begin prologue     sppose_z
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            calculate zero time gaussian wavepacket
!***
!***references
!***routines called
!***end prologue       sppose_z
  SUBROUTINE sppose_z(wave_function)
  IMPLICIT NONE
  LOGICAL                                        :: dollar
  INTEGER                                        :: ntot
  INTEGER, DIMENSION(3)                          :: n_val
  INTEGER                                        :: i, j, k, intkey
  COMPLEX*16, DIMENSION(n3d)                     :: wave_function
  ALLOCATE(c_loc(spdim))
  DO i=1,spdim
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
       IF (spdim >= 2) THEN
           n_val(2)=intkey(card,'number-of-y-superposed-states',1,' ')
           ntot=ntot*n_val(2)
           CALL intarr(card,'y-state-list',c_loc(2)%l,n_val(2),' ')
           CALL fparr(card,'y-state-coefficients',c_loc(2)%c,n_val(2),' ')
           CALL mk_phi(c_loc(2)%phi,grid(2)%eigvec_0,c_loc(2)%c,c_loc(2)%l, &
                       n_val(2),nphy(2),'y')
       END IF
       IF (spdim >= 3) THEN
           n_val(3)=intkey(card,'number-of-z-superposed-states',1,' ')
           ntot=ntot*n_val(3)
           CALL intarr(card,'z-state-list',c_loc(3)%l,n_val(3),' ')
           CALL fparr(card,'z-state-coefficients',c_loc(3)%c,n_val(3),' ')
           CALL mk_phi(c_loc(3)%phi,grid(3)%eigvec_0,c_loc(3)%c,c_loc(3)%l, &
                       n_val(3),nphy(3),'z')
      END IF
  END IF
  IF (spdim ==1 ) THEN
      CALL wfn_fil_1d_zd(wave_function,c_loc(1)%phi)
  ELSE IF( spdim == 2) THEN
      CALL wfn_fil_2d_zd(wave_function,c_loc(1)%phi,c_loc(2)%phi)
  ELSE IF( spdim == 3) THEN
      CALL wfn_fil_3d_zd(wave_function,c_loc(1)%phi,c_loc(2)%phi,c_loc(3)%phi)
  END IF
  DO i=1,spdim
     DEALLOCATE(c_loc(i)%c,c_loc(i)%phi,c_loc(i)%l)
  END DO
  DEALLOCATE(c_loc)
END SUBROUTINE sppose_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE initial_state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
