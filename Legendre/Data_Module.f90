!
MODULE Data_Module
!***begin prologue     Data_Module
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           time, finite element dvr, orthogonal polynomial
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            global shared variables for dvr library
!***description        this routine defines the global variables
!***                   and data needed for the dvrlib
!
!***references

!***routines called    
!***end prologue       Data_Module
  USE accuracy
  USE input_output
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!---------------------------------------------------------------------
!                   Some Constants and Units
!---------------------------------------------------------------------
!
   REAL(idp)                              ::          pi    = 3.141592653589793238462643383276D0 
   REAL(idp)                              ::          two_pi= 6.283185307179586476925286766552D0 
   REAL(idp)                              ::          zero    =  0.D0
   REAL(idp)                              ::          quarter = .25D0
   REAL(idp)                              ::          half    = .50D0
   REAL(idp)                              ::          third   = .33333333333333333333333333333D0
   REAL(idp)                              ::          fourth  = .25000000000000000000000000000D0
   REAL(idp)                              ::          fifth   = .20000000000000000000000000000D0
   REAL(idp)                              ::          sixth   = .16666666666666666666666666666D0
   REAL(idp)                              ::          seventh = .14285714285714285714285714285D0
   REAL(idp)                              ::          eighth  = .12500000000000000000000000000D0
   REAL(idp)                              ::          ninth   = .11111111111111111111111111111D0
   REAL(idp)                              ::          tenth   = .10000000000000000000000000000D0
   REAL(idp)                              ::          one     = 1.0D0
   REAL(idp)                              ::          two     = 2.0D0
   REAL(idp)                              ::          three   = 3.0D0
   REAL(idp)                              ::          four    = 4.0D0
   REAL(idp)                              ::          five    = 5.0D0
   REAL(idp)                              ::          six     = 6.0D0
   REAL(idp)                              ::          seven   = 7.0D0
   REAL(idp)                              ::          eight   = 8.0D0
   REAL(idp)                              ::          nine    = 9.0D0
   REAL(idp)                              ::          ten     = 10.D0
   REAL(idp)                              ::          nrzero  = 1.D-06
   INTEGER                                ::          int_zero      = 0
   INTEGER                                ::          int_one       = 1
   INTEGER                                ::          int_two       = 2
   INTEGER                                ::          int_three     = 3
   INTEGER                                ::          int_four      = 4
   INTEGER                                ::          int_five      = 5
   INTEGER                                ::          int_six       = 6
   INTEGER                                ::          int_seven     = 7
   INTEGER                                ::          int_eight     = 8
   INTEGER                                ::          int_nine      = 9
   INTEGER                                ::          int_ten       = 10
   INTEGER                                ::          int_eleven    = 11
   INTEGER                                ::          int_twelve    = 12
   INTEGER                                ::          int_thirteen  = 13
   INTEGER                                ::          int_fourteen  = 14
   INTEGER                                ::          int_fifteen   = 15
   INTEGER                                ::          int_sixteen   = 16
   INTEGER                                ::          int_seventeen = 17
   INTEGER                                ::          int_eighteen  = 18
   INTEGER                                ::          int_nineteen  = 19
   INTEGER                                ::          int_twenty    = 20
!
!                    hbar in joule-sec
!
   REAL(idp)                              ::          hbar = 1.054571596D-34
!
!                    electron mass in kg
!
   REAL(idp)                              ::          massau = 9.10938188D-31
!
!                    bohr radius in meters
!
   REAL(idp)                              ::          lenau = 5.291772083D-11
!
!                    time for an electron to make one bohr orbit in seconds
!
   REAL(idp)                              ::          timau    = 2.418884326D-17
   REAL(idp)                              ::          efieldau = 5.14220624D+11
   REAL(idp)                              ::          electric_field_to_intensity = 3.509338D+16
   REAL(idp)                              ::          peak_electric_field = .2849540283D-03
   REAL(idp)                              ::          pmass    = 1.67262158D-27
   REAL(idp)                              ::          massn2p  = 1.00137841887D0
   REAL(idp)                              ::          au_in_ev = 27.211396132D0
!
!
END MODULE Data_Module
