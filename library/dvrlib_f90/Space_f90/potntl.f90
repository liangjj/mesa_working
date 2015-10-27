!deck potntl
!***begin prologue     potntl
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            potentials
!***
!***references
!***routines called
!***end prologue       potntl
  SUBROUTINE potntl(v_coul,v_ang,pt,n)
  USE dvr_global,    ONLY    : inp, iout
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: v_coul, v_ang, pt
  CHARACTER(LEN=80)                      :: title      
  v_coul = 1.d0 / pt
  v_ang  = 1.d0/(pt*pt)
  if(prn(6)) then
     title='inverse r potential'
     call prntrm(title,v_coul,n,1,n,1,iout)
     title='inverse r**2 potential'
     call prntrm(title,v_ang,n,1,n,1,iout)
  end if
END SUBROUTINE potntl



