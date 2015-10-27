!deck Exponential_Time_Propagation
!***begin prologue     Exponential_Time_Propagation
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            
!***
!***references
!***routines called    iosys, util and mdutil
!***modules used
!***end prologue       Exponential_Time_Propagation
!
!
  Subroutine Exponential_Time_Propagation
  USE Time_Integrator_Derived_Types
  USE Runge_Kutta
  IMPLICIT NONE
  CHARACTER (LEN=1)                        :: itoc
  CHARACTER (LEN=80)                       :: chrkey
  LOGICAL                                  :: dollar
  LOGICAL                                  :: logkey
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
END Subroutine Exponential_Time_Propagation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
