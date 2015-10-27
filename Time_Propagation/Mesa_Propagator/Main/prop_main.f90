!deck prop_main
!**begin prologue     prop_main
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            main propagator program
!**description        calls routines to input data, construct FEDVR propagators
!**                   and then propagate the initial wave function.
!**                   
!**references
!**routines called      Name                    Location    
!                       ----                    --------
!                     Input_Prop             Prop_Sources
!                     Space_Prop             Prop_Sources
!                     Regional_Matrices      Prop_Modules:regional_module
!                     v_couple               Prop_Sources
!                     Arnoldi_Propagation    Prop_Modules:arnoldi_module
!                     SO_Sector_Time_Prop    Prop_Modules:split_operator_module
!                     SO_Propagation         Prop_Modules:split_operator_module
!
!**modules used         Name                    Location               
!                       ----                    --------
!                    dvrprop_global          Modules   
!                    arnoldi_module          Prop_Modules
!                    split_operator_module   Prop_Modules
!
!**end prologue       prop_main
  PROGRAM prop_main
  USE dvrprop_global
  USE arnoldi_module
  USE split_operator_module
  IMPLICIT NONE
  REAL*4                                   :: secnds
  REAL*4, DIMENSION(10)                    :: del_t
  INTEGER                                  :: iostat 
  CHARACTER (LEN=4096)                     :: ops
!
  CALL Drum
  CALL IOsys('read character options from rwf',-1,0,0,ops)
! Read in the Input
!
  time(1)=secnds(0.0)
  CALL Input_Prop
  time(2)=secnds(0.0)
  del_t(1) = time(2) - time(1)
!
! Calculate the Global FEDVR Matrices
!
  CALL Space_Prop
  time(3)=secnds(0.0)
  del_t(2) = time(3) - time(2)
!
! Compute the Regional Matrices from the FEDVR Global Matrices
!
  CALL Regional_Matrices
  time(4)=secnds(0.0)
  del_t(3) = time(4)-time(3)
!
! Get potential parameters.  The actual potential itself
! is calculated in the propagation step.
!
  CALL v_couple
  time(5)=secnds(0.0)
  del_t(4) = time(5)-time(4)
  WRITE(iout,1) del_t(1:4)
  WRITE(iout,4)
!
! ALLOCATE the Memory for the Matrices used to Store the Wavefunction
! and for the Potential.  Both the Arnoldi and Split Operator need these
! arrays.
!
! Here is where we have to begin to parallelize the code.
! Arrays have to be allocated to different processors and I am not sure of
! the best way to proceed.
!
  OPEN(UNIT=99,FILE='initial-wavefunction',ACCESS='sequential',            &
       FORM='unformatted',IOSTAT=IOSTAT,STATUS='unknown')
  IF(algorithm == 'arnoldi') THEN
!                                                                        
!                                                                         
!    Read in any data for the arnoldi iterations.                       
!                                                                         
     CALL Arnoldi_Data
!        
!    The Arnoldi method needs lots of memory in order to be efficient.
!    The big arrays are n3d*maxvec in size where maxvec is the maximum
!    number of Arnoldi vectors allowed per timestep.  The Split Operator
!    method is much less demanding and this is a big plus for many applications.
!      
     CALL Arnoldi_Driver
     time(6) = secnds(0.0)
     del_t(5) = time(6) - time(5)
     WRITE(iout,2) del_t(5)
     WRITE(iout,4)
  ELSE IF(algorithm == 'split_operator') THEN
     CALL Split_Operator_Driver
     time(6) = secnds(0.0)
     del_t(5) = time(6)-time(5)
     WRITE(iout,3) del_t(5)
     WRITE(iout,4)
  ELSE
    CALL lnkerr('algorithm error')
  END IF
  write(iout,*) 'Calling chainx'
  CALL chainx(0)
  write(iout,*) 'Called chainx'
1 FORMAT('***********************************************'                 &
            '*************************'                                    &
            /,10X,'Time to input basic data               = ',f15.8,/,     &
            /,10X,'Time to compute global FEDVR spatial ',                 &
            /,10x,'matrices                               = ',f15.8,/,     &
            /,10x,'Time to compute the regional spatial '                  &
            /,10x,'matrices                               = ',f15.8,/,     &
            /,10x,'Time to enter potential parameters     = ',f15.8)
2 FORMAT('***********************************************'                 &
            '*************************'                                    &
            /,10x,'Time for the Arnoldi propagation       = ',f15.8)
3 FORMAT('***********************************************'                 &
            '*************************'                                    &
            /,10x,'Time for the Split Operator Propagation = ',f15.8)
4 FORMAT('***********************************************'                 &
         '*************************')
  stop
END PROGRAM prop_main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
