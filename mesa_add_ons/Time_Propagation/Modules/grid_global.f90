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
   REAL*8         :: pi=3.141592653589793238462643D+00 
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
   REAL*8         ::           pmass    = 1.67262158D-27
   REAL*8         ::           massn2p  = 1.00137841887D0

!
!
  CHARACTER (LEN=80)                    :: keywrd, units, typpot
  CHARACTER (LEN=80)                    :: typwt, typarg, dentyp
  CHARACTER (LEN=80)                    :: typint, refwt, cpass, drctv
  CHARACTER (LEN=1600)                  :: card
  LOGICAL                               :: reuse, nodiag
  REAL*8,  DIMENSION(2)                 :: endpts
  INTEGER                               :: angmom
  REAL*8                                :: charge, mass
!
!
END MODULE grid_global
