!
MODULE Coulomb_Shared
!***begin prologue     Coulomb_Shared
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           coulomb wavefunctions
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            list of global shared variables for single center coulomb code
!***description        
!***                   
!
!***references

!***routines called    
!***end prologue       Coulomb_Shared
  USE accuracy
  USE input_output
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!---------------------------------------------------------------------
!                    shared variables used in lots of routines
!---------------------------------------------------------------------
  REAL(idp)                                          :: charge
  REAL(idp)                                          :: rmin
  REAL(idp)                                          :: rmax
  REAL(idp)                                          :: rdel
  REAL(idp)                                          :: energy
  REAL(idp)                                          :: rswtch
  REAL(idp)                                          :: eta
  REAL(idp)                                          :: twoeta
  REAL(idp)                                          :: zero = 0.d0
  REAL(idp)                                          :: one = 1.d0
  REAL(idp)                                          :: two = 2.d0
  REAL(idp)                                          :: sqrt2 = 1.414213562373095048801688724209d0
  COMPLEX(idp)                                       :: eye=(0.d0,1.d0)
  LOGICAl                                            :: drven
  LOGICAl                                            :: tstspl
  LOGICAl                                            :: nospln
  LOGICAl                                            :: noireg
!
!
!---------------------------------------------------------------------
!                    Keywords for Print Options
!---------------------------------------------------------------------
!
  CHARACTER (LEN=80), DIMENSION(12)           :: prnkey
  CHARACTER (LEN=80)                          :: secprn
  DATA prnkey /'sector_points', 'sector_factors',                         &
               'sector_polynomials', 'sector_matrices',                   &
               'global_points', 'global_polynomials',                     &
               'potential','global_matrices',                             &
               'hamiltonian', 'eigenvalues','eigenvectors', 'all' /
  DATA secprn / 'none' /
!
!
 CHARACTER (LEN=80),DIMENSION(12)             :: prloc
 LOGICAL, DIMENSION(12)                       :: prn
!
!
!***********************************************************************
!***********************************************************************
END MODULE Coulomb_Shared
