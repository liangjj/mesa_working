!************************************************************************************
!deck read_fourier
!***begin prologue     read_fourier
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       read_fourier
  SUBROUTINE Read_Fourier
  USE dvr_global
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                          :: ntest
  nreg=1
  nfix=0
  fix(1)=.false.
  fix(2)=.false.
  drop(1)=.false.
  drop(2)=.false.
  bcl=1
  bcr=1
!  edge(1)=-pi
!  edge(2)=pi
  edge(1)=0.d0
  edge(2)=2.d0*pi
  CALL fparr(card,'region_boundaries',edge,nreg+1,' ')
  CALL intarr(card,'number_of_points',n,nreg,' ')
!
!             if n is not odd fix it.
!
  ntest = n(1) - 2 * (n(1)/2)
  IF(ntest == 0 ) THEN
     n(1) = n(1) + 1
  END IF 
  npt(1)=n(1)
  nrq(1)=npt(1) 
  END SUBROUTINE Read_Fourier
