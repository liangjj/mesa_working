!deck Set_FEDVR_Defaults
!***begin prologue     Set_FEDVR_Defaults
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           Input, DVR
!***author             schneider, b. i.(nsf)
!***source
!***purpose            input for DVR
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
  Subroutine Set_FEDVR_Defaults
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  USE potential
  IMPLICIT NONE
  CHARACTER (LEN=1)                        :: itoc
  CHARACTER (LEN=80)                       :: chrkey
  LOGICAL                                  :: dollar
  LOGICAL                                  :: logkey
  LOGICAL                                  :: compute_regional_matrices
  INTEGER                                  :: intkey
  INTEGER                                  :: i

!
  IF ( dollar('$dvr_basis',card,cpass,inp) ) then
       spdim = intkey(card,'number_of_space_variables',1,' ')
       system='cartesian'
       units='atomic_units'
       hbar=1.d0
       mass=1.d0
       diag=.true.
       typke='dvr'
       proj=.false.
       diag_mod='none'
       keywrd='h0'
  END IF
  DO i = 1, spdim
     coord(i)=chrkey(card,'space_variable_'//itoc(i),'x',' ')
  END DO
END Subroutine Set_FEDVR_Defaults
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
