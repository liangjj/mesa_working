!
MODULE grid_global
!***begin prologue     grid_global
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           time, finite element dvr, orthogonal polynomial
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            global variables for dvr library
!***description        this routine defines the global variables
!***                   and data needed for the dvrlib
!
!***references

!***routines called    
!***end prologue       grid_global
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                   Some Constants and Units
!
   REAL*8         :: pi    = 3.141592653589793238462643383276D+00 
   REAL*8         :: two_pi= 6.283185307179586476925286766552D+00 
   REAL*8         :: zero    =  0.d0
   REAL*8         :: quarter = .25d0
   REAL*8         :: half    = .50d0
   REAL*8         :: one     = 1.0d0
   REAL*8         :: two     = 2.0d0
   REAL*8         :: three   = 3.0d0
   REAL*8         :: four    = 4.0d0
   REAL*8         :: five    = 5.0d0
   REAL*8         :: six     = 6.0d0
   REAL*8         :: seven   = 7.0d0
   REAL*8         :: eight   = 8.0d0
   REAL*8         :: nine    = 9.0d0
   REAL*8         :: ten     = 10.d0
   INTEGER        :: int_one    = 1
   INTEGER        :: int_two    = 2
   INTEGER        :: int_three  = 3
   INTEGER        :: int_four   = 4
   INTEGER        :: int_five   = 5
   INTEGER        :: int_six    = 6
   INTEGER        :: int_seven  = 7
   INTEGER        :: int_eight  = 8
   INTEGER        :: int_nine   = 9
   INTEGER        :: int_ten    = 10
   REAL*8         :: nrzero  = 1.d-06
!
!                    hbar in joule-sec
!
   REAL*8         ::           hbar = 1.054571596D-34
!
!                    electron mass in kg
!
   REAL*8         ::           massau = 9.10938188D-31
!
!                    bohr radius in meters
!
   REAL*8         ::           lenau = 5.291772083D-11
!
!                    time for an electron to make one bohr orbit in seconds
!
   REAL*8         ::           timau    = 2.418884326D-17
   REAL*8         ::           efieldau = 5.14220624D+11
   REAL*8         ::           electric_field_to_intensity = 3.509338D+16
   REAL*8         ::           peak_electric_field = .2849540283D-03
   REAL*8         ::           pmass    = 1.67262158D-27
   REAL*8         ::           massn2p  = 1.00137841887D0
   REAL*8         ::           au_in_ev = 27.211396132D0

!
!
  CHARACTER (LEN=80)                    :: keywrd
  CHARACTER (LEN=80)                    :: units
  CHARACTER (LEN=80)                    :: typpot
  CHARACTER (LEN=80)                    :: typwt
  CHARACTER (LEN=80)                    :: typarg
  CHARACTER (LEN=80)                    :: dentyp
  CHARACTER (LEN=80)                    :: typint
  CHARACTER (LEN=80)                    :: refwt
  CHARACTER (LEN=80)                    :: cpass
  CHARACTER (LEN=80)                    :: drctv
  CHARACTER (LEN=2400)                  :: card
  LOGICAL                               :: reuse
  LOGICAL                               :: nodiag
  LOGICAL                               :: reuse_sector_information
  LOGICAL                               :: unit_weight
  REAL*8,  DIMENSION(2)                 :: endpts
  INTEGER                               :: angmom
  REAL*8                                :: charge
  REAL*8                                :: mass
!
!
END MODULE grid_global
