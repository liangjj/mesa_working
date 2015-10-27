!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              MODULE initial_state_module
                              USE prop_global
                              USE dvrprop_global
                              USE dvr_shared
                              USE dvr_global
                              USE Matrix_Module
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
  INTEGER                                :: intkey, state
!
!       Read the data
!
  IF( dollar('$initial_state',card,title,inp) ) THEN
      i0stat=chrkey(card,'driver','unit_vector',' ')
      prnton=logkey(card,'print=on',.false.,' ')
      WRITE(iout,1) i0stat
  END IF
1 FORMAT(/,1X,'initial state data',                          &
         /,1X,'type initial state = ',a24)
END SUBROUTINE initial_state_data
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
  SUBROUTINE initial_vector_d(wave_function)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                   :: wave_function
  INTEGER                                :: IOSTAT
  INTEGER                                :: i
  CHARACTER (LEN=8)                      :: kewrd
  IF(i0stat == 'unit_vector' ) THEN
     wave_function = 1.d0
  ELSE IF(i0stat == 'random_vector' ) THEN
     call random_number(wave_function) 
  ELSE IF(i0stat.eq.'from_disk') THEN
     OPEN(UNIT=20,FILE='solution',ACCESS='sequential', FORM='unformatted', &
     IOSTAT=IOSTAT,STATUS='old')
     READ(20) wave_function(:)
     CLOSE(20)  
  ELSE IF(i0stat == 'right_hand_side') THEN
     wave_function = rhs_d/matrix_diagonal_d
  ELSE
     CALL lnkerr('error in initial state')
  END IF
  IF(prnton) THEN
     title='initial state'
     CALL prntfm(title,wave_function,n3d,1,n3d,1,iout,'e')
  END IF
  WRITE(99) wave_function
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
  SUBROUTINE initial_vector_z(wave_function)
  USE dvrprop_global,                    local_scratch => v_tot
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)               :: wave_function
  COMPLEX*16                             :: seed
  INTEGER                                :: i
  INTEGER                                :: IOSTAT
  REAL*8                                 :: t_ran
  CHARACTER (LEN=8)                      :: kewrd
  IF(i0stat == 'unit_vector' ) THEN
     wave_function(:) = (1.d0,0.d0)
  ELSE IF(i0stat == 'random_vector' ) THEN
     DO i=1,n3d
        call random_number(t_ran) 
        seed=cmplx(t_ran,t_ran)
        wave_function(i) = seed
     END DO
  ELSE IF(i0stat.eq.'from_disk') THEN
     OPEN(UNIT=20,FILE='eigenstates',ACCESS='sequential', FORM='unformatted', &
     IOSTAT=IOSTAT,STATUS='old')
     READ(20) local_scratch(:)
     CLOSE(20)      
     wave_function = local_scratch
  ELSE IF(i0stat == 'right_hand_side') THEN
     wave_function = rhs_z/matrix_diagonal_z
  ELSE
     CALL lnkerr('error in initial state')
  END IF
  IF(prnton) THEN
     title='initial state'
     CALL prntcm(title,wave_function,n3d,1,n3d,1,iout,'e')
  END IF
  WRITE(99) wave_function
END SUBROUTINE initial_vector_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE initial_state_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
