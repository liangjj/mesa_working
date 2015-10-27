!***********************************************************************
!deck Atomic_Main
!***begin prologue     Atomic_Main
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input, atomic
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for computing dvr basis sets,
!***                   kinetic energy, one-electron and two electron
!***                   radial integrals.
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       Atomic_Main
!***********************************************************************
  Program Atomic_Main
!***********************************************************************
  USE dvr_shared
  USE dvr_global
  USE dvrprop_global
  USE Atomic_Module
  IMPLICIT NONE
!
!     read the input by stopping at the keyword in the input stream.
  CALL Atomic_Input
  CALL Atomic_Basis
  CALL One_Electron_Integrals
  CALL Two_Electron_Integrals
!***********************************************************************
  END Program Atomic_Main
!***********************************************************************
