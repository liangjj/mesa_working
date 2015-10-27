!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     MODULE finite_element_propagator
!**begin prologue     finite_element_propagator
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            module containing the 1, 2 and 3D split operator
!**                   and arnoldi routines.  these routines drive the
!**                   propagation.
!**references
!**routines called
!**end prologue       finite_element_propagator
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     INTERFACE so_prop
             MODULE PROCEDURE so_prop_1d, so_prop_2d, so_prop_3d
                 END INTERFACE so_prop
                     INTERFACE arn_prop
             MODULE PROCEDURE arn_prop_1d, arn_prop_2d, arn_prop_3d
                 END INTERFACE arn_prop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck so_prop_1d
!**begin prologue     so_prop_1d
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
!**end prologue       so_prop_1d
  SUBROUTINE so_prop_1d(wave_function,scratch_vector)
  USE arnoldi_global_it
  USE dvrprop_global_it
  USE dvr_global
  USE dvr_shared
  USE normalize
  USE moment
  USE auto_correlation
  USE plot_wavefunction
  USE initial_state
  USE real_space_propagator
  USE psi_h_psi
  USE h_on_vector
  USE eigenvalues
  USE non_linear_potential
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1))              :: wave_function
  REAL*8, DIMENSION(nphy(1))              :: scratch_vector
  REAL*8                                  :: e_cur, eig_tst
  CHARACTER (LEN=2)                       :: itoc
  INTEGER                                 :: i, t, tst
  INTEGER                                 :: length, len
  REAL*8                                  :: norm
  COMPLEX*16                              :: auto
  REAL*8                                  :: total_t
!
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
      CALL cp_psi(wave_function,t,tim_pts(t))
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
      CALL rsp(wave_function,scratch_vector)  
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
         CALL eigen_solve(scratch_vector,wave_function,e_cur)
     ELSE
         call check_norm(wave_function,norm)
         wave_function = wave_function/sqrt(norm) 
         call h_v_d(wave_function,scratch_vector,1)
         call check_energy(wave_function,scratch_vector,e_cur)
     END IF
     eig_tst=abs(e_cur-eig_old)
     IF( eig_tst <= con ) THEN
!
!        End the propagation
!
         CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
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
!
!
END SUBROUTINE so_prop_1d
!deck so_prop_2d
!**begin prologue     so_prop_2d
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
!**end prologue       so_prop_2d
  SUBROUTINE so_prop_2d(wave_function,scratch_vector)
  USE arnoldi_global_it
  USE dvrprop_global_it
  USE dvr_global
  USE dvr_shared
  USE normalize
  USE moment
  USE auto_correlation
  USE plot_wavefunction
  USE initial_state
  USE real_space_propagator
  USE eigenvalues
  USE psi_h_psi
  USE h_on_vector
  USE non_linear_potential
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1))         :: wave_function
  REAL*8, DIMENSION(nphy(2),nphy(1))         :: scratch_vector
  CHARACTER (LEN=2)                          :: itoc
  INTEGER                                    :: i, t, tst
  INTEGER                                    :: length, len
  REAL*8                                     :: norm, rtemp
  REAL*8                                     :: e_cur, eig_tst
  COMPLEX*16                                 :: auto
  REAL*8                                     :: total_t
!
  eig_old=0.d0
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
      CALL cp_psi(wave_function,t,tim_pts(t))
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
      CALL rsp(wave_function,scratch_vector)    
      tst = t - plot_step * ( t /plot_step )
!
!  
!        compare the approximate and exact solutions where possible
!
     IF( tst == 0 ) THEN
!
!       Form the total solution and compute the autocorrelation
!       function.
        CALL plot_psi(wave_function,t)
     END IF
     IF (e_method == 'exponential') then
         CALL eigen_solve(scratch_vector,                             &
                          wave_function,f_1,e_cur)
     ELSE
         call check_norm(wave_function,norm)
         wave_function = wave_function/sqrt(norm) 
         call h_v_d(wave_function,scratch_vector,1)
         call check_energy(wave_function,                             &
                           scratch_vector,f_1,e_cur)
     END IF
     eig_tst=abs(e_cur-eig_old)
     IF( eig_tst <= con ) THEN
!
!        End the propagation
!
         CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
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
!
!
END SUBROUTINE so_prop_2d
!deck so_prop_3d
!**begin prologue     so_prop_3d
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
!**end prologue       so_prop_3d
  SUBROUTINE so_prop_3d(wave_function,scratch_vector)
  USE arnoldi_global_it
  USE dvrprop_global_it
  USE dvr_global
  USE dvr_shared
  USE normalize
  USE moment
  USE auto_correlation
  USE plot_wavefunction
  USE initial_state
  USE real_space_propagator
  USE eigenvalues
  USE psi_h_psi
  USE h_on_vector
  USE non_linear_potential
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))       :: wave_function
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))       :: scratch_vector
  CHARACTER (LEN=2)                                :: itoc
  INTEGER                                          :: i, t, tst
  INTEGER                                          :: iostat, length, len
  REAL*8                                           :: norm, rtemp
  REAL*8                                           :: e_cur, eig_tst
  COMPLEX*16                                       :: auto
  REAL*8                                           :: t_int, total_t
!
  eig_old=0.d0
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
      CALL cp_psi(wave_function,t,tim_pts(t))
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
      CALL rsp(wave_function,scratch_vector) 
      tst = t - plot_step * ( t /plot_step )
!
!  
!        compare the approximate and exact solutions where possible
!
     IF( tst == 0 ) THEN
!
!       Form the total solution and compute the autocorrelation
!       function.
        CALL plot_psi(wave_function,t)
     END IF
     IF (e_method == 'exponential') then
         CALL eigen_solve(scratch_vector,                                  &
                          wave_function,f_1,f_2,e_cur)
     ELSE
         call check_norm(wave_function,norm)
         wave_function = wave_function/sqrt(norm) 
         call h_v_d(wave_function,scratch_vector,1)
         call check_energy(wave_function,                                  &
                           scratch_vector,f_1,f_2,e_cur)
     END IF
     eig_tst=abs(e_cur-eig_old)
     IF( eig_tst <= con ) THEN
!
!        End the propagation
!
         CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
         RETURN
     ELSE
         eig_old=e_cur
     END IF  
     call iosys('write real solution to bec',n3d,                          &
                 wave_function,0,' ')
  END DO
!
!     End the propagation
!
  CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
  call check_norm(wave_function,norm)
  wave_function = wave_function/sqrt(norm) 
  CALL plot_psi(wave_function,t-1) 
!
!
END SUBROUTINE so_prop_3d
!
!deck arn_prop_1d
!**begin prologue     arn_prop_1d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Arnoldi, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Arnoldi_Prop_f90
!**purpose            time dependent schroedinger equation
!**                   using Arnoldi method with finite difference or
!**                   dvr space representation.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       arn_prop_1d_real
  SUBROUTINE arn_prop_1d(wave_function,scratch_vector)
  USE arnoldi_driver
  USE arnoldi_global_it
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  USE normalize
  USE moment
  USE auto_correlation
  USE plot_wavefunction
  USE initial_state
  USE real_space_propagator
  USE psi_h_psi
  USE h_on_vector
  USE eigenvalues
  USE non_linear_potential
  USE pack_and_read_h
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1))              :: wave_function
  REAL*8, DIMENSION(nphy(1))              :: scratch_vector
  CHARACTER (LEN=2)                       :: itoc
  INTEGER                                 :: i, t, tst
  INTEGER                                 :: bigs, bigv
  INTEGER                                 :: length, len
  REAL*8                                  :: t_int, total_t
  REAL*8                                  :: norm
  ALLOCATE(buf(1))
  ALLOCATE( buf(1)%d(nphy(1)),            &
            buf(1)%hbuf(nphy(1)*nphy(1)), &
            buf(1)%hibuf(2,nphy(1)*nphy(1)))
  CALL pack_h(1)
  eig_old=0.d0
  DO  t=1,ntreg
      keywrd=itoc(t)
      LEN=length(keywrd)
      keywrd='t'//keywrd(1:LEN)
      LEN=length(keywrd)
      keywrd='$v0('//keywrd(1:LEN)//')'
      t0=tim_pts(t)
      t1=tim_pts(t+1)
      deltat=t1-t0
  
!        Calculate the time dependent perturbation.
!        It consists of a space and a time part.

      v_tot = 0.d0
      CALL v_tim
      CALL pert
!  
!     Initialize the wavefunction at $t_{0}$ or 
!     read in its value from the disk at $t_{i-1}$.
!
      CALL cp_psi(wave_function,t,tim_pts(t))
!  
!     Calculate the non-linear potential if present.
!  
      IF(vtyp(2) /= 'none') THEN
         CALL v_nl(wave_function)
      END IF
!
!        Propagate to the next time using the short time formula,
!\begin{equation}
!    \Xi(x,y,z,t_{i-1} + {\delta}t) = \big [ exp( - i H(x,y,z,t_{i-1}) {\delta}t )
!                                -1 \big ] \Psi_{0}(x,y,z,t_{i-1})
!                   \nonumber
!\end{equation}
!        This is done in the Arnoldi code by performing an effective reduction
!        of the Hamiltonian to diagonal form using the value of the wavefunction
!        from the previous step as an initial guess.
!
     CALL arnoldi_d(wave_function,scratch_vector)
!  
!        compare the approximate and exact solutions where possible
     tst = t - plot_step * ( t /plot_step )  
     IF(tst == 0) THEN
!
!       Form the total solution and compute the autocorrelation
!       function.
!
        CALL plot_psi(wave_function,t)
     END IF
     call iosys('write real solution to bec',n3d,                   &
                 wave_function,0,' ')
! Form the total solution and compute the autocorrelation
! function.
  
     call check_norm(wave_function,norm)
     IF( cntrl == 'done' ) THEN
!
!        End the propagation
!
         CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
         call check_norm(wave_function,norm)
         wave_function = wave_function/sqrt(norm)
         CALL plot_psi(wave_function,t-1) 
         DEALLOCATE( buf(1)%d,buf(1)%hbuf,buf(1)%hibuf)
         DEALLOCATE(buf)
         RETURN
     END IF  
  END DO
!
END SUBROUTINE arn_prop_1d
!deck arn_prop_2d
!**begin prologue     arn_prop_2d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Arnoldi, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Arnoldi_Prop_f90
!**purpose            time dependent schroedinger equation
!**                   using Arnoldi method with finite difference or
!**                   dvr space representation.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       arn_prop_2d
  SUBROUTINE arn_prop_2d(wave_function,scratch_vector)
  USE arnoldi_driver
  USE arnoldi_global_it
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  USE normalize
  USE moment
  USE auto_correlation
  USE plot_wavefunction
  USE initial_state
  USE real_space_propagator
  USE psi_h_psi
  USE h_on_vector
  USE eigenvalues
  USE non_linear_potential
  USE pack_and_read_h
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1))        :: wave_function
  REAL*8, DIMENSION(nphy(2),nphy(1))        :: scratch_vector
  CHARACTER (LEN=2)                         :: itoc
  INTEGER                                   :: i, t, tst
  INTEGER                                   :: bigs, bigv
  INTEGER                                   :: length, len
  REAL*8                                    :: t_int, total_t
  REAL*8                                    :: norm
  ALLOCATE(buf(2))
  ALLOCATE( buf(1)%d(nphy(1)),                                   &
            buf(1)%hbuf(nphy(1)*nphy(1)),                        &
            buf(1)%hibuf(2,nphy(1)*nphy(1)),                     &
            buf(2)%d(nphy(2)),                                   &
            buf(2)%hbuf(nphy(2)*nphy(2)),                        &
            buf(2)%hibuf(2,nphy(2)*nphy(2)))
  CALL pack_h(1)
  CALL pack_h(2)
  eig_old=0.d0
  DO  t=1,ntreg
      keywrd=itoc(t)
      LEN=length(keywrd)
      keywrd='t'//keywrd(1:LEN)
      LEN=length(keywrd)
      keywrd='$v0('//keywrd(1:LEN)//')'
      t0=tim_pts(t)
      t1=tim_pts(t+1)
      deltat=t1-t0
  
!        Calculate the time dependent perturbation.
!        It consists of a space and a time part.

      v_tot = 0.d0
      CALL v_tim
      CALL pert
!  
!     Initialize the wavefunction at $t_{0}$ or 
!     read in its value from the disk at $t_{i-1}$.
!
      CALL cp_psi(wave_function,t,tim_pts(t))
!  
!     Calculate the non-linear potential if present.
!  
      IF(vtyp(2) /= 'none') THEN
         CALL v_nl(wave_function)
      END IF
!
!        Propagate to the next time using the short time formula,
!\begin{equation}
!    \Xi(x,y,z,t_{i-1} + {\delta}t) = \big [ exp( - i H(x,y,z,t_{i-1}) {\delta}t )
!                                -1 \big ] \Psi_{0}(x,y,z,t_{i-1})
!                   \nonumber
!\end{equation}
!        This is done in the Arnoldi code by performing an effective reduction
!        of the Hamiltonian to diagonal form using the value of the wavefunction
!        from the previous step as an initial guess.
!
      CALL arnoldi_d(wave_function,scratch_vector)
!  
!        compare the approximate and exact solutions where possible
     tst = t - plot_step * ( t /plot_step )  
     IF(tst == 0) THEN
!
!       Form the total solution and compute the autocorrelation
!       function.
!
        CALL plot_psi(wave_function,t)
     END IF
! Form the total solution and compute the autocorrelation
! function.
     call iosys('write real solution to bec',n3d,                &
                 wave_function,0,' ')  
     call check_norm(wave_function,norm)
     IF( cntrl == 'done' ) THEN
!
!        End the propagation
!
         CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
         call check_norm(wave_function,norm)
         wave_function = wave_function/sqrt(norm)
         CALL plot_psi(wave_function,t-1) 
         DEALLOCATE( buf(1)%d,buf(1)%hbuf,buf(1)%hibuf)
         DEALLOCATE(buf)
         RETURN
     END IF  
   END DO
!
!
END SUBROUTINE arn_prop_2d
!deck arn_prop_3d
!**begin prologue     arn_prop_3d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Arnoldi, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Arnoldi_Prop_f90
!**purpose            time dependent schroedinger equation
!**                   using Arnoldi method with finite difference or
!**                   dvr space representation.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       arn_prop_3d
  SUBROUTINE arn_prop_3d(wave_function,scratch_vector)
  USE arnoldi_driver
  USE arnoldi_global_it
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  USE normalize
  USE moment
  USE auto_correlation
  USE plot_wavefunction
  USE initial_state
  USE real_space_propagator
  USE psi_h_psi
  USE h_on_vector
  USE eigenvalues
  USE non_linear_potential
  USE pack_and_read_h
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))            &
                                                 :: wave_function
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))            &
                                                 :: scratch_vector
  CHARACTER (LEN=2)                              :: itoc
  INTEGER                                        :: i, t, tst
  INTEGER                                        :: bigs, bigv
  INTEGER                                        :: length, len
  REAL*8                                         :: t_int, total_t
  REAL*8                                         :: norm
  ALLOCATE(buf(3))
  ALLOCATE( buf(1)%d(nphy(1)),                                   &
            buf(1)%hbuf(nphy(1)*nphy(1)),                        &
            buf(1)%hibuf(2,nphy(1)*nphy(1)),                     &
            buf(2)%d(nphy(2)),                                   &
            buf(2)%hbuf(nphy(2)*nphy(2)),                        &
            buf(2)%hibuf(2,nphy(2)*nphy(2)),                     &
            buf(3)%d(nphy(3)),                                   &
            buf(3)%hbuf(nphy(3)*nphy(3)),                        &
            buf(3)%hibuf(2,nphy(3)*nphy(3)))
  CALL pack_h(1)
  CALL pack_h(2)
  CALL pack_h(3)
  eig_old=0.d0
  DO  t=1,ntreg
      keywrd=itoc(t)
      LEN=length(keywrd)
      keywrd='t'//keywrd(1:LEN)
      LEN=length(keywrd)
      keywrd='$v0('//keywrd(1:LEN)//')'
      t0=tim_pts(t)
      t1=tim_pts(t+1)
      deltat=t1-t0
  
!        Calculate the time dependent perturbation.
!        It consists of a space and a time part.

      v_tot = 0.d0
      CALL v_tim
      CALL pert
!  
!     Initialize the wavefunction at $t_{0}$ or 
!     read in its value from the disk at $t_{i-1}$.
!
      CALL cp_psi(wave_function,t,tim_pts(t))
!  
!     Calculate the non-linear potential if present.
!  
      IF(vtyp(2) /= 'none') THEN
         CALL v_nl(wave_function)
      END IF
!
!        Propagate to the next time using the short time formula,
!\begin{equation}
!    \Xi(x,y,z,t_{i-1} + {\delta}t) = \big [ exp( - i H(x,y,z,t_{i-1}) {\delta}t )
!                                -1 \big ] \Psi_{0}(x,y,z,t_{i-1})
!                   \nonumber
!\end{equation}
!        This is done in the Arnoldi code by performing an effective reduction
!        of the Hamiltonian to diagonal form using the value of the wavefunction
!        from the previous step as an initial guess.
!
      CALL arnoldi_d(wave_function,scratch_vector)
!  
!        compare the approximate and exact solutions where possible
     tst = t - plot_step * ( t /plot_step )  
     IF(tst == 0) THEN
!
!       Form the total solution and compute the autocorrelation
!       function.
!
        CALL plot_psi(wave_function,t)
     END IF
! Form the total solution and compute the autocorrelation
! function.
     call iosys('write real solution to bec',n3d,                &
                 wave_function,0,' ')  
     call check_norm(wave_function,norm)
     IF( cntrl == 'done' ) THEN
!
!        End the propagation
!
         CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
         call check_norm(wave_function,norm)
         wave_function = wave_function/sqrt(norm)
         CALL plot_psi(wave_function,t-1) 
         DEALLOCATE( buf(1)%d,buf(1)%hbuf,buf(1)%hibuf)
         DEALLOCATE(buf)
         RETURN
     END IF  
   END DO
!
!
END SUBROUTINE arn_prop_3d
END MODULE finite_element_propagator
