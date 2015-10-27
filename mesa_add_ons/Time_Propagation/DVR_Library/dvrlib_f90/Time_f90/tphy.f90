!deck tphy.f
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            physical hamiltonian.
!***
!***references
!***routines called
!***end prologue       tphy

  SUBROUTINE tphy(tin,twtin,ptin,dptin,ddptin,hin,  &
                  t0out,tout,twtout,ptout,dptout,ddptout, &
                  hout,m,n,prn)
  USE dvr_global,    ONLY  : iout
  IMPLICIT NONE
  INTEGER                                :: m, n
  REAL*8, DIMENSION(m)                   :: tin, twtin
  REAL*8, DIMENSION(m,m)                 :: ptin, dptin, ddptin, hin
  REAL*8                                 :: t0out
  REAL*8, DIMENSION(n)                   :: tout, twtout
  REAL*8, DIMENSION(n,n)                 :: ptout, dptout, ddptout, hout
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  t0out=tin(1)
  CALL copy(tin(2),tout,n)
  CALL copy(twtin(2),twtout,n)
  CALL mmove(ptin(2,2),ptout,n,n,m,n)
  CALL mmove(dptin(2,2),dptout,n,n,m,n)
  CALL mmove(ddptin(2,2),ddptout,n,n,m,n)
  CALL mmove(hin(2,2),hout,n,n,m,n)
  IF(prn) THEN
     title='initial grid'
     CALL prntrm(title,tin,m,1,m,1,iout)
     title='initial weights'
     CALL prntrm(title,twtin,m,1,m,1,iout)
     title='initial functions'
     CALL prntrm(title,ptin,m,m,m,m,iout)
     title='initial first derivatives'
     CALL prntrm(title,dptin,m,m,m,m,iout)
     title='initial second derivatives'
     CALL prntrm(title,ddptin,m,m,m,m,iout)
     title='initial hamiltonian'
     CALL prntrm(title,hin,m,m,m,m,iout)
     title='final grid'
     CALL prntrm(title,tout,n,1,n,1,iout)
     title='final weights'
     CALL prntrm(title,twtout,n,1,n,1,iout)
     title='final functions'
     CALL prntrm(title,ptout,n,n,n,n,iout)
     title='final first derivatives'
     CALL prntrm(title,dptout,n,n,n,n,iout)
     title='final second derivatives'
     CALL prntrm(title,ddptout,n,n,n,n,iout)
     title='final hamiltonian'
     CALL prntrm(title,hout,n,n,n,n,iout)
   END IF
   RETURN
END SUBROUTINE tphy




