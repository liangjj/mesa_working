!************************************************************************************
!deck read_legendre
!***begin prologue     read_legendre
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       read_legendre
  SUBROUTINE Read_Legendre
  USE dvr_global
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                               :: intkey
  LOGICAL                               :: logkey
  nreg=1
  m_val=intkey(card,'legendre-m',0,' ')
  write(iout,1) m_val
!
!     This will do what needs to be done for one region
!     If there is more than one region the user needs to call
!     the grid_parameters subroutine and explicitly set things up.
!
  skip=logkey(card,'read-grid-parameters',.false.,' ')
  IF (skip) THEN
      CALL Read_Grid_Parameters
  ELSE
      CALL intarr(card,'number-of-points',n,nreg,' ')
      npt(1)=n(1)
      nrq(1)=npt(1) 
      IF(m_val == 0 ) THEN
         nfix=0
         fix(1)=.false.
         fix(2)=.false.
         drop(1)=.false.
         drop(2)=.false.
         bcl=1
         bcr=1
         edge(1)=-1.d0
         edge(2)=1.d0
      ELSE 
         nfix=2
         fix(1)=.true.
         fix(2)=.true.
         drop(1)=.true.
         drop(2)=.true.
         bcl=0
         bcr=0
      END IF 
  END IF
1 Format(/,1x,'legendre m value = ',i2)

!*************************************************************************************
  END SUBROUTINE Read_Legendre
