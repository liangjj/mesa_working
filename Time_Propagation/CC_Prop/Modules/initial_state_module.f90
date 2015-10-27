!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              MODULE initial_state_module
                              USE prop_global
                              USE dvrprop_global
                              USE dvr_shared
                              USE dvr_global
                              USE plot_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**begin prologue     initial_state_module
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate the spatial part of the initial wavepacket.
!**references
!**                   
!**routines called    c_vect, rad_paket, gauss_paket, sppose, z_proj
!**end prologue       initial_state_module
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              INTERFACE gauss_packet
              MODULE PROCEDURE gauss_packet_d,                          &
                               gauss_packet_z
                          END INTERFACE gauss_packet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              INTERFACE z_proj
              MODULE PROCEDURE z_proj_d,                                &
                               z_proj_z                       
                          END INTERFACE z_proj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck initial_state_data
!**begin prologue     initial_state_data
!**date written       040706   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            get data to compute initial state.
!**references
!**                   
!**routines called    c_vect, rad_paket, gauss_paket, sppose, z_proj
!**end prologue       initial_state_data
!
  SUBROUTINE initial_state_data
  IMPLICIT NONE
  LOGICAL                                :: dollar, logkey
  CHARACTER (LEN=2)                      :: itoc
  CHARACTER (LEN=80)                     :: chrkey
  INTEGER                                :: intkey
!
!       Read the data
!
  IF( dollar('$initial_state',card,title,inp) ) THEN
      i0stat=chrkey(card,'initial_state','unperturbed_state_vector',' ')
      prnton=logkey(card,'print=on',.false.,' ')
      WRITE(iout,1) i0stat
  END IF
  IF(i0stat == 'gaussian_pulse') THEN
     CALL gaussian_data
  ELSE IF (i0stat.eq.'unperturbed_state_vector') THEN
     state=intkey(card,'initial_state',0,' ')
     write(iout,2) state
  ELSE IF(i0stat.eq.'perturbed_state_vector') THEN
     state=intkey(card,'initial_state',0,' ')
     write(iout,2) state
  ELSE IF(i0stat.eq.'unit_vector') THEN
     write(iout,3)     
  ELSE IF(i0stat.eq.'random_vector') THEN
     write(iout,4)      
  ELSE IF(i0stat.eq.'from_standard_disk'.or.i0stat.eq.'from_mesa_disk') THEN
     write(iout,5)      
     state=intkey(card,'initial_state',0,' ')
  ELSE
     CALL lnkerr('error in initial state')
  END IF
1 FORMAT(/,1X,'initial state data',                             &
         /,1X,'type initial state = ',a24)
2 FORMAT(/,5x,'the initial state = ',i3)
3 FORMAT(/,5x,'the initial state = unit_vector')
4 FORMAT(/,5x,'the initial state = random_vector')
5 FORMAT(/,5x,'the initial state = from_disk')
END SUBROUTINE initial_state_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck gaussian_data
!***begin prologue     gaussian_data
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            input data for gaussian calculation.
!                      
!***
!***references
!***routines called
!***end prologue       gaussian_data
!
  SUBROUTINE gaussian_data
  IMPLICIT NONE                          
  INTEGER                                :: intkey
  REAL*8                                 :: fpkey
  LOGICAL                                :: dollar
  CHARACTER(LEN=80)                      :: chrkey
!
!
  IF ( dollar('$initial_state',card,title,inp) ) THEN
       CALL fparr(card,'alpha',alpha,spdim,' ')
       CALL fparr(card,'sigma',sigma,spdim,' ')
       CALL fparr(card,'x_0',x_0,spdim,' ')
       CALL fparr(card,'beta',beta,spdim,' ')
       powr=intkey(card,'power',0,' ')
       typ_pak=chrkey(card,'type_radial_packet', &
                           'exponential',' ')
  ELSE
       write(iout,1)
       call lnkerr
  END IF
  sigma = sigma/sqrt(alpha)
  alpha = 1.d0
1 format(/,1x,'This is a test for a separable free wavepacket', &
         /,1x,'Current calculation does not qualify')
END SUBROUTINE gaussian_data
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
!**routines called    c_vect, rad_paket, gauss_paket, sppose, z_proj
!**end prologue       initial_state
!
  SUBROUTINE initial_vector_d(wave_function,file_directory)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                   :: wave_function
  INTEGER                                :: i , nwfn
  INTEGER                                :: len
  INTEGER                                :: lenth
  REAL*8                                 :: ddot, norm, nfac
  CHARACTER(LEN=8)                       :: kewrd
  CHARACTER(LEN=*)                       :: file_directory
  SAVE norm
  nwfn=size(wave_function,1)
!
!     The allowable right hand sides are; one, a state vector of the
!     unperturbed Hamiltonian, a specified superposition of 
!     states of the unperturbed Hamiltonian or a Gaussian 
!     pulse with or without a velocity component.  The routine could 
!     easily be modified to do more.
  IF(i0stat.eq.'unperturbed_state_vector') THEN
     IF( spdim == 1) THEN
         energy=grid(1)%eigv_0(state+1)
         wave_function(:) = grid(1)%eigvec_0(:,state+1)   
     ELSE IF (spdim == 2 ) THEN
         CALL c_vect_2d_d(wave_function,etemp,itemp,state+1)
     ELSE IF (spdim == 3) THEN
         CALL c_vect_3d_d(wave_function,etemp,itemp,state+1)
     END IF
  ELSE IF(i0stat.eq.'perturbed_state_vector') THEN
     IF( spdim == 1) THEN
         energy=grid(1)%eigv(state+1)
         wave_function(:) = grid(1)%eigvec(:,state+1)   
     END IF 
  ELSE IF(i0stat == 'gaussian_pulse') THEN
     CALL nr_packet(norm)
     CALL gauss_packet(wave_function,norm)
  ELSE IF(i0stat == 'superpose' ) THEN
     CALL sppose(wave_function)
  ELSE IF(i0stat == 'unit_vector' ) THEN
     wave_function = 1.d0
  ELSE IF(i0stat == 'random_vector' ) THEN
     call random_number(wave_function) 
  ELSE IF(i0stat.eq.'from_standard_disk') THEN
     Write(iout,*) 'Reading Initial Vector From Standard Disk'
     REWIND(20)
     READ(20) kewrd
     DO i=0,state
        READ(20) wave_function(:)
     END DO
  ELSE IF(i0stat.eq.'from_mesa_disk') THEN
     len=lenth(file_directory)
     Write(iout,*) 'Reading Initial Vector From Time_Dependent_Data'
     Call IOsys('open time_dependent_data as unknown',0,0,0,                    &
                 file_directory(1:len)//'/'//'Time_Dependent_Data')
     Call IOsys('read real ground_state_eigenvector from time_dependent_data',  &
                 nwfn,wave_function,0,' ')
     Call IOsys('close time_dependent_data',0,0,0,' ')
  ELSE
     CALL lnkerr('error in initial state')
  END IF
  IF(i0stat /='from_standard_disk'.or.i0stat /='from_mesa_disk') THEN  
      norm=ddot(nwfn,wave_function,1,wave_function,1)
      write(iout,2) norm
      wave_function = wave_function/sqrt(norm) 
  END IF
  IF(prnton) THEN
     title='initial state'
     CALL print_psi(wave_function)
  END IF
  WRITE(99) wave_function
1 FORMAT(/,1X,'initial state construction at first time = ',e15.8,  &
         /,1X,'driver                                   = ',a24)
2 FORMAT(/,5x,'Normalization Integral = ',e15.8)
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
!**routines called    c_vect, rad_paket, gauss_paket, sppose, z_proj
!**end prologue       initial_state
!
  SUBROUTINE initial_vector_z(wave_function,file_directory)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)               :: wave_function
  COMPLEX*16                             :: cdotc
  COMPLEX*16                             :: seed
  INTEGER                                :: i , nwfn
  INTEGER                                :: len
  INTEGER                                :: lenth
  REAL*8                                 :: norm, nfac
  REAL*8                                 :: t_ran
  CHARACTER(LEN=8)                       :: kewrd
  CHARACTER(LEN=*)                       :: file_directory
  REAL*8, DIMENSION(:), ALLOCATABLE      :: r_temp
  SAVE norm
  nwfn=size(wave_function,1)
!    The allowable right hand sides are; one, a state vector of the
!    unperturbed Hamiltonian, a specified superposition of 
!    states of the unperturbed Hamiltonian or a Gaussian 
!    pulse with a velocity component.  The routine could easily be 
!    modified to do more.
  IF(i0stat.eq.'unperturbed_state_vector') then
     IF( spdim == 1) THEN
         energy=grid(1)%eigv_0(state+1)
         wave_function(:) = grid(1)%eigvec_0(:,state+1)   
     ELSE IF (spdim == 2 ) THEN
         CALL c_vect_2d_z(wave_function,etemp,itemp,state+1)
     ELSE IF (spdim == 3) THEN
         CALL c_vect_3d_z(wave_function,etemp,itemp,state+1)
     END IF
  ELSE IF(i0stat.eq.'perturbed_state_vector') THEN
     IF( spdim == 1) THEN
         energy=grid(1)%eigv(state+1)
         wave_function(:) = grid(1)%eigvec(:,state+1)   
     END IF 
  ELSE IF(i0stat == 'gaussian_pulse') THEN
     CALL nr_packet(norm)
     CALL gauss_packet(wave_function,norm)
  ELSE IF(i0stat == 'superpose' ) THEN
     CALL sppose(wave_function)
  ELSE IF(i0stat == 'unit_vector' ) THEN
     wave_function(:) = (1.d0,0.d0)
  ELSE IF(i0stat == 'random_vector' ) THEN
     DO i=1,nwfn
        call random_number(t_ran) 
        seed=cmplx(t_ran,t_ran)
        wave_function(i) = seed
     END DO
  ELSE IF(i0stat.eq.'from_standard_disk') THEN
     Write(iout,*) 'Reading Initial Vector From Standard Disk'
     REWIND(20)
     READ(20) kewrd
     IF(kewrd == 'real' ) THEN
        ALLOCATE(r_temp(nwfn))
        DO i=0,state
           READ(20) r_temp(:)  
        END DO
        DO i=1,nwfn
           wave_function(i) = r_temp(i)
        END DO
        DEALLOCATE(r_temp)
     ELSE IF(kewrd == 'complex' ) THEN
        DO i=0,state
           READ(20) wave_function  
        END DO
     END IF
  ELSE IF(i0stat.eq.'from_mesa_disk') THEN
     Write(iout,*) 'Reading Initial Vector From Time_Dependent_Data'
     len=lenth(file_directory)
     Call IOsys('open time_dependent_data as unknown',0,0,0,                    &
                 file_directory(1:len)//'/'//'Time_Dependent_Data')
     Call IOsys('read real ground_state_eigenvector from time_dependent_data',  &
                 2*nwfn,wave_function,0,' ')
     Call IOsys('close time_dependent_data',0,0,0,' ')
  ELSE
     CALL lnkerr('error in initial state')
  END IF
  IF(prnton) THEN
     title='initial state'
     CALL print_psi(wave_function)
  END IF
  WRITE(99) wave_function
1 FORMAT(/,1X,'initial state construction at first time = ',e15.8,  &
         /,1X,'driver                                   = ',a24)
2 FORMAT(/,5x,'Normalization Integral = ',e15.8)
END SUBROUTINE initial_vector_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**begin prologue     gaussian_wave_packet
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate zero time gaussian wavepacket
!**references
!**routines called
!**end prologue       gaussian_wave_packet
! \center Gaussian Wave Packet
! \par
!\begin{eqnarray}
! \Psi({\bf r}) & = N *
!               & { ( x -x _0 ) }^{n}_x exp( - {\alpha}_x \frac{ ( x -x_0 )^2 }
!                                                           {2 { {\sigma}_x }^2 } )
!                                                            & \nonumber \\
!               & { ( y - y_0 ) }^{n}_y exp( - {\alpha}_y \frac{ ( y -y_0 )^2 } 
!                                                           {2 { {\sigma}_y }^2 } )
!                                                            & \nonumber \\
!               & { (z -z _0 ) }^{n}_z exp( - {\alpha}_z \frac{ ( z -z_0 )^2 }
!                                                           {2 { {\sigma}_z }^2 } )
!                                                            &
! \end{eqnarray}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE gauss_packet_d(wave_function,norm)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)           :: wave_function
  REAL*8                         :: norm
  INTEGER                        :: i
  ALLOCATE(zloc(spdim))
  DO i=1, spdim
     ALLOCATE(zloc(i)%z_d(nphy(i)))
  END DO
  WRITE(iout,1)
  wave_function(:) = 0.d0
  IF ( spdim == 1 ) THEN
       WRITE(iout,2)
       WRITE(iout,3)
       WRITE(iout,4) alpha(1)
       WRITE(iout,5) sigma(1)
       WRITE(iout,6) x_0(1)
       WRITE(iout,7) powr(1)
       WRITE(iout,8) beta(1)
  ELSE IF (spdim == 2) THEN
       WRITE(iout,9)
       WRITE(iout,3)
       WRITE(iout,10) (alpha(i),i=1,2)
       WRITE(iout,11) (sigma(i),i=1,2)
       WRITE(iout,12) (x_0(i),i=1,2)
       WRITE(iout,13) (powr(i),i=1,2)
       WRITE(iout,14) (beta(i),i=1,2)
  ELSE IF (spdim == 3) THEN
       WRITE(iout,15)
       WRITE(iout,3)
       WRITE(iout,16) (alpha(i),i=1,2)
       WRITE(iout,17) (sigma(i),i=1,2)
       WRITE(iout,18) (x_0(i),i=1,2)
       WRITE(iout,19) (powr(i),i=1,2)
       WRITE(iout,20) (beta(i),i=1,2)
  END IF
  DO i=1,spdim
     CALL z_proj(grid(i)%pt,grid(i)%wt,                                    &
                 alpha(i),sigma(i),x_0(i),powr(i),beta(i),                 &
                 zloc(i)%z_d,nphy(i))
  END DO
  IF (spdim == 1 ) THEN
      CALL wfn_fil_1d_dd(wave_function,zloc(1)%z_d)
  ELSE IF( spdim == 2) THEN
      CALL wfn_fil_2d_dd(wave_function,zloc(1)%z_d,zloc(2)%z_d)
  ELSE IF( spdim == 3) THEN
      CALL wfn_fil_3d_dd(wave_function,zloc(1)%z_d,zloc(2)%z_d,zloc(3)%z_d)
  END IF
  wave_function(:) = norm * wave_function(:)
  DO i=1,spdim
     DEALLOCATE(zloc(i)%z_d)
  END DO
  DEALLOCATE(zloc)
1  FORMAT(/,5X,'initial wavepacket at t=0')
2  FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,    &
    'psi = exp[ - alpha_x * ( x - x_0 )**2 / (2*(sigma_x)**2)'   &
              ' - i*beta_x*(x-x_0) ]' )
3 FORMAT(/,1X,'gaussian wave packet parameters')
4 FORMAT(/,1X,'alpha_x    = ',e15.8)
5 FORMAT(/,1X,'sigma_x    = ',e15.8)
6 FORMAT(/,1X,'x_0        = ',e15.8)
7 FORMAT(/,1X,'n_x        = ',i2)
8 FORMAT(/,1X,'beta_x     = ',e15.8)
9 FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,    &
    'psi = exp[ - alpha_x * ( x - x_0 )**2 / (2*(sigma_x)**2)'   &
              ' - i*beta_x*(x-x_0) ]',/,5X,                      &
    '      exp[ - alpha_y * ( y - y_0 )**2 / (2*(sigma_y)**2)'   &
              ' - i*beta_y*(y-y_0) ]' )
10 FORMAT(/,1X,'alpha_x    = ',e15.8,1X,'alpha_y   = ',e15.8)
11 FORMAT(/,1X,'sigma_x    = ',e15.8,1X,'sigma_y   = ',e15.8)
12 FORMAT(/,1X,'x_0        = ',e15.8,1X,'y_0       = ',e15.8)
13 FORMAT(/,1X,'n_x        = ',i2,   1X,'n_y       = ',i2 )
14 FORMAT(/,1X,'beta_x     = ',e15.8,1X,'beta_y    = ',e15.8 )
15 FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,    &
    'psi = exp[ - alpha_x * ( x - x_0 )**2 / (2*(sigma_x)**2)'   &
              ' - i*beta_x*(x-x_0) ]',/,5X,                      &
    '      exp[ - alpha_y * ( y - y_0 )**2 / (2*(sigma_y)**2)'   &
              ' - i*beta_y*(y-y_0) ]',/,5x,                      &
    '      exp[ - alpha_z * ( z - z_0 )**2 / (2*(sigma_z)**2)'   &
              ' - i*beta_z*(z-z_0) ]' )
16 FORMAT(/,1X,'alpha_x    = ',e15.8,1X,'alpha_y   = ',e15.8,1X,'alpha_z   = ',e15.8)
17 FORMAT(/,1X,'sigma_x    = ',e15.8,1X,'sigma_y   = ',e15.8,1X,'sigma_z   = ',e15.8)
18 FORMAT(/,1X,'x_0        = ',e15.8,1X,'y_0       = ',e15.8,1X,'z_0       = ',e15.8)
19 FORMAT(/,1X,'n_x        = ',i2,   1X,'n_y       = ',i2,   1X,'n_y       = ',i2)
20 FORMAT(/,1X,'beta_x     = ',e15.8,1X,'beta_y    = ',e15.8,'beta_z       = ',e15.8)
  END SUBROUTINE gauss_packet_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE gauss_packet_z(wave_function,norm)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)       :: wave_function
  REAL*8                         :: norm
  INTEGER                        :: i
  ALLOCATE(zloc(spdim))
  DO i=1, spdim
     ALLOCATE(zloc(i)%z_z(nphy(i)))
  END DO
  WRITE(iout,1)
  wave_function(:) = 0.d0
  IF ( spdim == 1 ) THEN
       WRITE(iout,2)
       WRITE(iout,3)
       WRITE(iout,4) alpha(1)
       WRITE(iout,5) sigma(1)
       WRITE(iout,6) x_0(1)
       WRITE(iout,7) powr(1)
       WRITE(iout,8) beta(1)
  ELSE IF (spdim == 2) THEN
       WRITE(iout,9)
       WRITE(iout,3)
       WRITE(iout,10) (alpha(i),i=1,2)
       WRITE(iout,11) (sigma(i),i=1,2)
       WRITE(iout,12) (x_0(i),i=1,2)
       WRITE(iout,13) (powr(i),i=1,2)
       WRITE(iout,14) (beta(i),i=1,2)
  ELSE IF (spdim == 3) THEN
       WRITE(iout,15)
       WRITE(iout,3)
       WRITE(iout,16) (alpha(i),i=1,2)
       WRITE(iout,17) (sigma(i),i=1,2)
       WRITE(iout,18) (x_0(i),i=1,2)
       WRITE(iout,19) (powr(i),i=1,2)
       WRITE(iout,20) (beta(i),i=1,2)
  END IF
  DO i=1,spdim
     CALL z_proj(grid(i)%pt,grid(i)%wt,                                    &
                 alpha(i),sigma(i),x_0(i),powr(i),beta(i),                 &
                 zloc(i)%z_z,nphy(i))
  END DO
  IF (spdim == 1 ) THEN
      CALL wfn_fil_1d_zz(wave_function,zloc(1)%z_z)
  ELSE IF( spdim == 2) THEN
      CALL wfn_fil_2d_zz(wave_function,zloc(1)%z_z,zloc(2)%z_z)
  ELSE IF( spdim == 3) THEN
      CALL wfn_fil_3d_zz(wave_function,zloc(1)%z_z,zloc(2)%z_z,zloc(3)%z_z)
  END IF
  wave_function(:) = norm * wave_function(:)
  DO i=1,spdim
     DEALLOCATE(zloc(i)%z_z)
  END DO
  DEALLOCATE(zloc)
1  FORMAT(/,5X,'initial wavepacket at t=0')
2  FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,    &
    'psi = exp[ - alpha_x * ( x - x_0 )**2 / (2*(sigma_x)**2)'   &
              ' - i*beta_x*(x-x_0) ]' )
3 FORMAT(/,1X,'gaussian wave packet parameters')
4 FORMAT(/,1X,'alpha_x    = ',e15.8)
5 FORMAT(/,1X,'sigma_x    = ',e15.8)
6 FORMAT(/,1X,'x_0        = ',e15.8)
7 FORMAT(/,1X,'n_x        = ',i2)
8 FORMAT(/,1X,'beta_x     = ',e15.8)
9 FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,    &
    'psi = exp[ - alpha_x * ( x - x_0 )**2 / (2*(sigma_x)**2)'   &
              ' - i*beta_x*(x-x_0) ]',/,5X,                      &
    '      exp[ - alpha_y * ( y - y_0 )**2 / (2*(sigma_y)**2)'   &
              ' - i*beta_y*(y-y_0) ]' )
10 FORMAT(/,1X,'alpha_x    = ',e15.8,1X,'alpha_y   = ',e15.8)
11 FORMAT(/,1X,'sigma_x    = ',e15.8,1X,'sigma_y   = ',e15.8)
12 FORMAT(/,1X,'x_0        = ',e15.8,1X,'y_0       = ',e15.8)
13 FORMAT(/,1X,'n_x        = ',i2,   1X,'n_y       = ',i2 )
14 FORMAT(/,1X,'beta_x     = ',e15.8,1X,'beta_y    = ',e15.8 )
15 FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,    &
    'psi = exp[ - alpha_x * ( x - x_0 )**2 / (2*(sigma_x)**2)'   &
              ' - i*beta_x*(x-x_0) ]',/,5X,                      &
    '      exp[ - alpha_y * ( y - y_0 )**2 / (2*(sigma_y)**2)'   &
              ' - i*beta_y*(y-y_0) ]',/,5x,                      &
    '      exp[ - alpha_z * ( z - z_0 )**2 / (2*(sigma_z)**2)'   &
              ' - i*beta_z*(z-z_0) ]' )
16 FORMAT(/,1X,'alpha_x    = ',e15.8,1X,'alpha_y   = ',e15.8,1X,'alpha_z   = ',e15.8)
17 FORMAT(/,1X,'sigma_x    = ',e15.8,1X,'sigma_y   = ',e15.8,1X,'sigma_z   = ',e15.8)
18 FORMAT(/,1X,'x_0        = ',e15.8,1X,'y_0       = ',e15.8,1X,'z_0       = ',e15.8)
19 FORMAT(/,1X,'n_x        = ',i2,   1X,'n_y       = ',i2,   1X,'n_y       = ',i2)
20 FORMAT(/,1X,'beta_x     = ',e15.8,1X,'beta_y    = ',e15.8,'beta_z       = ',e15.8)
  END SUBROUTINE gauss_packet_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wfn_fil_1d_dd(wave_function,z_1)
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1))            :: wave_function
  REAL*8, DIMENSION(nphy(1))            :: z_1
  wave_function(:) = z_1(:)
  END SUBROUTINE wfn_fil_1d_dd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wfn_fil_1d_zd(wave_function,z_1)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(1))        :: wave_function
  REAL*8, DIMENSION(nphy(1))            :: z_1
  wave_function(:) = z_1(:)
  END SUBROUTINE wfn_fil_1d_zd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wfn_fil_1d_zz(wave_function,z_1)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(1))            :: wave_function
  COMPLEX*16, DIMENSION(nphy(1))            :: z_1
  wave_function(:) = z_1(:)
  END SUBROUTINE wfn_fil_1d_zz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wfn_fil_2d_dd(wave_function,z_1,z_2)
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1))            :: wave_function
  REAL*8, DIMENSION(nphy(1))                    :: z_1
  REAL*8, DIMENSION(nphy(2))                    :: z_2
  INTEGER                                       :: i, j
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        wave_function(j,i) = z_1(i) * z_2(j)
     END DO
  END DO
  END SUBROUTINE wfn_fil_2d_dd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wfn_fil_2d_zd(wave_function,z_1,z_2)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(2),nphy(1))        :: wave_function
  REAL*8, DIMENSION(nphy(1))                    :: z_1
  REAL*8, DIMENSION(nphy(2))                    :: z_2
  INTEGER                                       :: i, j
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        wave_function(j,i) = z_1(i) * z_2(j)
     END DO
  END DO
  END SUBROUTINE wfn_fil_2d_zd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wfn_fil_2d_zz(wave_function,z_1,z_2)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(2),nphy(1))        :: wave_function
  COMPLEX*16, DIMENSION(nphy(1))                :: z_1
  COMPLEX*16, DIMENSION(nphy(2))                :: z_2
  INTEGER                                       :: i, j
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        wave_function(j,i) = z_1(i) * z_2(j)
     END DO
  END DO
  END SUBROUTINE wfn_fil_2d_zz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wfn_fil_3d_dd(wave_function,z_1,z_2,z_3)
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))    :: wave_function
  REAL*8, DIMENSION(nphy(1))                    :: z_1
  REAL*8, DIMENSION(nphy(2))                    :: z_2
  REAL*8, DIMENSION(nphy(3))                    :: z_3
  INTEGER                                       :: i, j, k
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        DO k=1,nphy(3)
           wave_function(k,j,i) = z_1(i) * z_2(j) * z_3(k)
        END DO
     END DO
  END DO
  END SUBROUTINE wfn_fil_3d_dd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wfn_fil_3d_zd(wave_function,z_1,z_2,z_3)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1))    :: wave_function
  REAL*8, DIMENSION(nphy(1))                        :: z_1
  REAL*8, DIMENSION(nphy(2))                        :: z_2
  REAL*8, DIMENSION(nphy(3))                        :: z_3
  INTEGER                                           :: i, j, k
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        DO k=1,nphy(3)
           wave_function(k,j,i) = z_1(i) * z_2(j) * z_3(k)
        END DO
     END DO
  END DO
  END SUBROUTINE wfn_fil_3d_zd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wfn_fil_3d_zz(wave_function,z_1,z_2,z_3)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1))  :: wave_function
  COMPLEX*16, DIMENSION(nphy(1))                  :: z_1
  COMPLEX*16, DIMENSION(nphy(2))                  :: z_2
  COMPLEX*16, DIMENSION(nphy(3))                  :: z_3
  INTEGER                                       :: i, j, k
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        DO k=1, nphy(3)
           wave_function(k,j,i) = z_1(i) * z_2(j) * z_3(k)
        END DO
     END DO
  END DO
  END SUBROUTINE wfn_fil_3d_zz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck z_proj_d
!***begin prologue     z_proj_d
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            tabulate gaussian wavepacket on the grid.
!***                   on grid
!***references
!***routines called
!***end prologue       z_proj_d
  SUBROUTINE z_proj_d(q,wt,alpha,sigma,x_0,n_p,beta,z,n)
  IMPLICIT NONE
  INTEGER                                :: n, n_p
  REAL*8, DIMENSION(n)                   :: q, wt
  REAL*8                                 :: alpha, sigma, x_0, beta
  REAL*8, DIMENSION(n)                   :: z
  INTEGER                                :: i
  IF(system == 'cartesian') THEN
     z = q**n_p * EXP( - alpha * ( q - x_0 ) * ( q - x_0 )  &
                        / ( 2.d0 * sigma * sigma))
  ELSE
     z = q**n_p * EXP( - alpha * ( q - x_0 ) * ( q - x_0 ) )
  END IF
  z = sqrt(wt) * z 
END SUBROUTINE z_proj_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck z_proj_z
!***begin prologue     z_proj_z
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            tabulate gaussian wavepacket
!***                   on grid
!***references
!***routines called
!***end prologue       z_proj_z
  SUBROUTINE z_proj_z(q,wt,alpha,sigma,x_0,n_p,beta,z,n)
  IMPLICIT NONE
  INTEGER                                :: n, n_p
  REAL*8, DIMENSION(n)                   :: q, wt
  REAL*8                                 :: alpha, sigma, x_0, beta
  COMPLEX*16, DIMENSION(n)               :: z
  INTEGER                                :: i
  IF(system == 'cartesian') THEN
     z = q**n_p * EXP( - alpha * ( q - x_0 ) * ( q - x_0 )  &
                        / ( 2.d0 * sigma * sigma))
     z = z * exp( - eye * beta * ( q - x_0) )
  ELSE
     z = q**n_p * EXP( - alpha * ( q - x_0 ) * ( q - x_0 ) )
     z = z * exp( - eye * beta * ( q - x_0) )
  END IF
  z = sqrt(wt) * z 
END SUBROUTINE z_proj_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck c_vect_2d_d.f
!**begin prologue     c_vect_2d_d
!**date written       040607   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           hamiltonian
!**author             schneider, barry (nsf)
!**source             2-dim
!**purpose            Sort 2-D eigenvectors according energy and place
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
  DO  ii=2,count
      i=ii-1
      k=i
      tmp=etemp(i)
      i1=itemp(i,1)
      j1=itemp(i,2)
      DO  j=ii,count
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
!**purpose            Sort 2-D eigenvectors according energy and place
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
  DO  ii=2,count
      i=ii-1
      k=i
      tmp=etemp(i)
      i1=itemp(i,1)
      j1=itemp(i,2)
      DO  j=ii,count
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
!**purpose            Sort 3-D eigenvectors according energy and place
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
  DO  ii=2,count
      i=ii-1
      k=i
      tmp=etemp(i)
      i1=itemp(i,1)
      j1=itemp(i,2)
      k1=itemp(i,3)
      DO  j=ii,count
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
!**purpose            Sort 3-D eigenvectors according energy and place
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
  DO  ii=2,count
      i=ii-1
      k=i
      tmp=etemp(i)
      i1=itemp(i,1)
      j1=itemp(i,2)
      k1=itemp(i,3)
      DO  j=ii,count
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
!***purpose            calculate a fixed superposition of eigenstates
!***                   as an initial wave packet.
!***references
!***routines called
!***end prologue       sppose_d
  SUBROUTINE sppose_d(wave_function)
  IMPLICIT NONE
  LOGICAL                                        :: dollar
  INTEGER                                        :: ntot
  INTEGER, DIMENSION(3)                          :: n_val
  INTEGER                                        :: i, j, k, intkey
  REAL*8, DIMENSION(:)                           :: wave_function
  wave_function = 0.d0
  IF ( dollar('$states',card,title,inp) ) THEN
       ntot=1
       n_val(1)=intkey(card,'number_of_x_superposed_states',1,' ')
       ntot=ntot*n_val(1)
       CALL intarr(card,'x_state_list',c_loc(1)%l,n_val(1),' ')
       CALL fparr(card,'x_state_coefficients',c_loc(1)%c,n_val(1),' ')
       CALL mk_phi(c_loc(1)%phi,grid(1)%eigvec_0,c_loc(1)%c,c_loc(1)%l, &
                   n_val(1),nphy(1),'x')
       IF (spdim >= 2) THEN
           n_val(2)=intkey(card,'number_of_y_superposed_states',1,' ')
           ntot=ntot*n_val(2)
           CALL intarr(card,'y_state_list',c_loc(2)%l,n_val(2),' ')
           CALL fparr(card,'y_state_coefficients',c_loc(2)%c,n_val(2),' ')
           CALL mk_phi(c_loc(2)%phi,grid(2)%eigvec_0,c_loc(2)%c,c_loc(2)%l, &
                       n_val(2),nphy(2),'y')
       END IF
       IF (spdim >= 3) THEN
           n_val(3)=intkey(card,'number_of_z_superposed_states',1,' ')
           ntot=ntot*n_val(3)
           CALL intarr(card,'z_state_list',c_loc(3)%l,n_val(3),' ')
           CALL fparr(card,'z_state_coefficients',c_loc(3)%c,n_val(3),' ')
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
!***purpose            calculate a fixed superposition of eigenstates
!***                   as the initial wave packet.
!***references
!***routines called
!***end prologue       sppose_z
  SUBROUTINE sppose_z(wave_function)
  IMPLICIT NONE
  LOGICAL                                        :: dollar
  INTEGER                                        :: ntot
  INTEGER, DIMENSION(3)                          :: n_val
  INTEGER                                        :: i, j, k, intkey
  COMPLEX*16, DIMENSION(:)                       :: wave_function
  wave_function = 0.d0
  IF ( dollar('$states',card,title,inp) ) THEN
       ntot=1
       n_val(1)=intkey(card,'number_of_x_superposed_states',1,' ')
       ntot=ntot*n_val(1)
       CALL intarr(card,'x_state_list',c_loc(1)%l,n_val(1),' ')
       CALL fparr(card,'x_state_coefficients',c_loc(1)%c,n_val(1),' ')
       CALL mk_phi(c_loc(1)%phi,grid(1)%eigvec_0,c_loc(1)%c,c_loc(1)%l, &
                   n_val(1),nphy(1),'x')
       IF (spdim >= 2) THEN
           n_val(2)=intkey(card,'number_of_y_superposed_states',1,' ')
           ntot=ntot*n_val(2)
           CALL intarr(card,'y_state_list',c_loc(2)%l,n_val(2),' ')
           CALL fparr(card,'y_state_coefficients',c_loc(2)%c,n_val(2),' ')
           CALL mk_phi(c_loc(2)%phi,grid(2)%eigvec_0,c_loc(2)%c,c_loc(2)%l, &
                       n_val(2),nphy(2),'y')
       END IF
       IF (spdim >= 3) THEN
           n_val(3)=intkey(card,'number_of_z_superposed_states',1,' ')
           ntot=ntot*n_val(3)
           CALL intarr(card,'z_state_list',c_loc(3)%l,n_val(3),' ')
           CALL fparr(card,'z_state_coefficients',c_loc(3)%c,n_val(3),' ')
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
END SUBROUTINE sppose_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE initial_state_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
