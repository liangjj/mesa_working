!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     MODULE so_propagator
                     USE arnoldi_global
                     USE dvrprop_global
                     USE dvr_global
                     USE dvr_shared
                     USE normalize
                     USE moment
                     USE auto_correlation
                     USE plot_wavefunction
                     USE spatial_wavefunction
                     USE initial_state
                     USE so_exponentiation
                     USE h_on_vector
                     USE non_linear_potential
                     USE pack_and_read_h
!**begin prologue     so_propagator
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            module containing the 1, 2 and 3D split operator
!**                   routines.  these routines drive the
!**                   propagation.
!**references
!**routines called
!**end prologue       finite_element_propagator
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     INTERFACE so_propagation
             MODULE PROCEDURE so_propagation_d,             &
                              so_propagation_z
                 END INTERFACE so_propagation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck so_propagation_d
!**begin prologue     so_propagation_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, real-space, split-operator, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             FEDVR_Prop
!**purpose            time dependent schroedinger equation
!**                   using split-operator method with 
!**                   finite difference or dvr space representation.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       so_propagation_d
  SUBROUTINE so_propagation_d(wave_function,scratch_vector)
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d)              :: wave_function
  REAL*8, DIMENSION(n3d)              :: scratch_vector
  CHARACTER (LEN=2)                   :: itoc
  INTEGER                             :: i, t, tst
  INTEGER                             :: length, len
  REAL*8                              :: norm
  COMPLEX*16                          :: auto
  REAL*8                              :: total_t
!
  IF(typke == 'dvr'.OR.typke == 'packed') THEN
     IF (e_method == 'hamiltonian') THEN
!
!           Calculate Energy in propagation using Hamiltonian
!
         ALLOCATE(buf(spdim))
         DO i=1,spdim
            ALLOCATE( buf(i)%d(nphy(i)),              &
                      buf(i)%hbuf(nphy(i)*nphy(i)),   &
                      buf(i)%hibuf(2,nphy(i)*nphy(i)))
            CALL pack_h(i)
         END DO
     END IF
  END IF
  WRITE(iout,1)
  WRITE(iout,2)
  WRITE(iout,1)
  eig_old=0.0d0
  DO  t=1,ntreg
      keywrd=itoc(t)
      LEN=length(keywrd)
      keywrd='t'//keywrd(1:LEN)
      LEN=length(keywrd)
      keywrd='$v0('//keywrd(1:LEN)//')'
      t0=tim_pts(t)
      t1=tim_pts(t+1)
  
!        Calculate the time dependent perturbation.  It consists of a
!        space and a time part.
!
      v_tot = 0.d0
      CALL v_tim
!
!        In subroutine pert, the one-body potentials are added to
!        any more complicated interaction.  Note that, the one-body
!        potentials may contain the diagonal contributions from the one-body
!        kinetic energy.
!
      CALL pert
!  
!     Initialize the wavefunction at $t_{0}$ or 
!     read in its value from the disk at $t_{i-1}$. 
!
      CALL initial_vector(wave_function,t,tim_pts(t))
!  
!     Calculate the non-linear potential if present.
!  
      IF(vtyp(2) /= 'none') THEN
         CALL v_nl(wave_function)
      END IF
!
!     Since the potential is, in general, time-dependent, the diagonal
!     propagator must be recomputed each time.
!
!     the main real space product driving routine
!
      CALL so_exp(wave_function,scratch_vector)  
!
      tst = t - plot_step * ( t /plot_step )
!
!  
!        compare the approximate and exact solutions where possible
!
     IF( tst == 0 ) THEN
!
!       Form the total solution and compute the autocorrelation
!       function.
!
        CALL plot_psi(wave_function,t)
     END IF
     IF (e_method == 'exponential') then
         CALL gs_energy(scratch_vector,wave_function,e_cur)
     ELSE
         call check_norm(wave_function,norm)
         wave_function = wave_function/sqrt(norm) 
         call h_v_d(wave_function,scratch_vector,1)
         call check_gs_energy(wave_function,scratch_vector,e_cur)
     END IF
     write(iout,3) t, norm, e_cur
     eig_tst=abs(e_cur-eig_old)
     IF( eig_tst <= con ) THEN
!
!        End the propagation
!
         CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
         WRITE(iout,1)
         WRITE(iout,4)
         WRITE(iout,1)
         RETURN
     ELSE
         eig_old=e_cur
     END IF  
     call iosys('write real solution to bec',n3d,                     &
                 wave_function,0,' ')
  END DO
!
!     End the propagation
!
  CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
  call check_norm(wave_function,norm)
  wave_function = wave_function/sqrt(norm)
  CALL plot_psi(wave_function,t-1) 
  IF(typke == 'dvr'.OR.typke == 'packed') THEN
     IF (e_method == 'hamiltonian') THEN
         DO i=1,spdim
            DEALLOCATE( buf(i)%d,buf(i)%hbuf,buf(i)%hibuf )
         END DO
         DEALLOCATE(buf)
     END IF
  END IF
!
1 FORMAT('***********************************************'                 &
         '*************************')
2 FORMAT(/,20x,'Begin Split Operator Propagation')
3 FORMAT(/,5x,'Time Step = ',i6,/,10x,'Normalization = ',e15.8,1x,     &
              'Energy    = ',e15.8)
4 FORMAT(/,20x,'End Split Operator Propagation')
!
END SUBROUTINE so_propagation_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck so_propagation_z
!**begin prologue     so_propagation_z
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, real-space, split-operator, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             FEDVR_Prop
!**purpose            time dependent schroedinger equation
!**                   using split-operator method with 
!**                   finite difference or dvr space representation.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       so_propagation_z
  SUBROUTINE so_propagation_z(wave_function,scratch_vector)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n3d)          :: wave_function
  COMPLEX*16, DIMENSION(n3d)          :: scratch_vector
  CHARACTER (LEN=2)                   :: itoc
  INTEGER                             :: i, t, tst
  INTEGER                             :: length, len, count
  REAL*8                              :: norm, rtemp
  COMPLEX*16                          :: auto
  REAL*8                              :: total_t
!
  WRITE(iout,1)
  WRITE(iout,2)
  WRITE(iout,1)
  count=0
  DO  t=1,ntreg
      keywrd=itoc(t)
      LEN=length(keywrd)
      keywrd='t'//keywrd(1:LEN)
      LEN=length(keywrd)
      keywrd='$v0('//keywrd(1:LEN)//')'
      t0=tim_pts(t)
      t1=tim_pts(t+1)
  
!        Calculate the time dependent perturbation.  It consists of a
!        space and a time part.
      v_tot = 0.d0
      CALL v_tim
!
!        In subroutine pert, the one-body potentials are added to
!        any more complicated interaction.  Note that, the one-body
!        potentials contain the diagonal contributions from the one-body
!        kinetic energy.
!
      CALL pert
!  
!     Initialize the wavefunction at $t_{0}$ or 
!     read in its value from the disk at $t_{i-1}$.  The first time through
!     the initial state is stored in soln which is needed later.
!
      CALL initial_vector(wave_function,t,tim_pts(t))
!  
!     Calculate the non-linear potential if present.
!  
      IF(vtyp(2) /= 'none') THEN
         CALL v_nl(wave_function)
      END IF
!
!     Since the potential is, in general, time-dependent, the diagonal
!     propagator must be recomputed each time.
!
!      The input initial state is psi.  The output is chi
!
      CALL so_exp(wave_function,scratch_vector)   
      call iosys('write real solution to bec',2*n3d,wave_function,0,' ')
!
!  
!        compare the approximate and exact solutions where possible
!
     IF( count == plot_step ) THEN
        count=0
        IF(i0stat == 'gaussian-pulse') THEN
           call calculate_moment(wave_function,tim_pts(t+1))
        ELSE IF(i0stat == 'perturbed-state-vector') THEN
!           CALL evolve(wave_function,scratch_vector,tim_pts(t+1))
        END IF
!       Form the total solution and compute the autocorrelation
!       function.
        CALL spatial_psi(wave_function,scratch_vector)
        CALL plot_psi(scratch_vector,t)
        CALL auto_correlation_function(auto,scratch_vector,     &
                                       wave_function)
        rtemp = 1.d0 - auto * conjg(auto)
        call check_norm(wave_function,norm)
        write(iout,1) t, norm, rtemp
     END IF
     count=count+1
  END DO
  WRITE(iout,1)
  WRITE(iout,4)
  WRITE(iout,1)
!
!     End the propagation
!
  CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
!
!
1 FORMAT('***********************************************'              &
         '*************************')
2 FORMAT(/,20x,'Begin Split Operator Propagation')
3 FORMAT(/,5x,'Time Step = ',i6,/,10x,                                  &
         /,5x,'Survival Probablity = ', e15.8)
4 FORMAT(/,20x,'End Split Operator Propagation')
END SUBROUTINE so_propagation_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE so_propagator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
