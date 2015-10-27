!deck Lanczos_Input.f
!***begin prologue     Lanczos_Input
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           input, propagation
!***author             schneider, b. i.(nsf)
!***source
!***purpose            input for close-coupled time propagation
!***
!***references
!***routines called    iosys, util and mdutil
!***modules used
!**                    Name                 Location
!                      ----                 -------
!
!                   dvrprop_global          Modules
!                   dvr_shared              Modules
!                   dvr_global              Modules
!
!***end prologue       input
!
!
  Subroutine Lanczos_Input(matrix_source)
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  USE Lanczos_Module
  IMPLICIT NONE
  CHARACTER (LEN=1)                        :: itoc
  CHARACTER (LEN=80)                       :: chrkey
  CHARACTER (LEN=8)                        :: matrix_source
  LOGICAL                                  :: dollar, logkey
  LOGICAL                                  :: propagation_test
  REAL*8                                   :: fpkey
  INTEGER                                  :: input, output
  INTEGER                                  :: intkey, i, j
  COMMON /io/ input, output
!
!  Get the input and output file numbers which appear in the 
!  Module input_output.f90 and put them into input and output which appears
!  in the ONLY common block in the code.  This is needed to pass into library
!  routines and for no other purpose.
!
  input = inp
  output = iout
!
!  Open the input and output files
!
  OPEN(input,file='Lanczos.inp',status='old')
  OPEN(output,file='Lanczos.out',status='unknown')
  WRITE(iout,*)
!
!  Look for the keyword $Lanczos_Setup in the input file and read all variables
!  that are associated with this keyword.
!
  spdim=1
  units='atomic-units'
  hbar=1.d0
  mass=1.d0
  IF ( dollar('$lanczos_setup',card,cpass,inp) ) then
!
       matrix_source=chrkey(card,'matrix_source','internal',' ')
!
!       Print options are set here.  Taking the default(none) will give minimal
!       print which is ok most of the time but not when debugging or getting
!       a feel for how the code runs.  Using the all-details option gives lots
!       of information.
!
       pr_main(1)=chrkey(card,'print','none',' ')
       log_main(1)=.false.
       IF(pr_main(1)=='all_details') THEN
          log_main(1)=.true.
       END IF
!
!
!
  ELSE
       write(iout,2)
       stop
  END IF  
!
!
!
1 FORMAT(/,20X,'Solve Lanczos Problem')
2 FORMAT(/,5x,'no basis card section')
END Subroutine Lanczos_Input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
