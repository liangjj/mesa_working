!deck rhofun
!***begin prologue     rhofun
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            model densities
!***
!***references
!***routines called
!***end prologue       rhofun
  SUBROUTINE rhofun(rho,pt,n)
  USE dvr_global,    ONLY    : inp, iout, dentyp
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: rho, pt
  CHARACTER(LEN=80)                      :: title      
  if (dentyp == 'exponential') then
      rho = pt*exp(-pt)
  else if (dentyp == 'harmonic') then
      rho = pt*exp(-pt*pt)
  else if (dentyp == 'one') then
      rho = pt
  else if (dentyp == 'linear') then
      rho = pt*pt
  end if
  if(prn(6)) then
     title='density'
     call prntrm(title,rho,n,1,n,1,iout)
  end if
END SUBROUTINE rhofun



