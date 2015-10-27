!deck tlobato.f
!***begin prologue     tlobato
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            To get the DVR points and weights
!                      in the time variable.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       tlobato

  SUBROUTINE tlobato(pt,wt,f,df,ddf,h,nphy)
  USE dvr_prnt
  USE dvr_global,   ONLY : typwt, n, edge 
  IMPLICIT NONE
  INTEGER nphy
  INTEGER hin, tin,twtin, ptin, dptin, ddptin, need
  REAL*8, dimension(nphy)        :: wt
  REAL*8, dimension(nphy+1)      :: pt
  REAl*8, dimension(nphy,nphy)   :: f, df, ddf, h
!
! Scratch array
!
  REAL*8, dimension(:), allocatable      :: scr
!
!    Define the local scratch which will be partitioned out and then destroyed.
!    I have done this as a single allocated array and then passed the indices.
  hin=1
  tin=hin+n(1)*n(1)
  twtin=tin+n(1)
  ptin=twtin+n(1)
  dptin=ptin+n(1)*n(1)
  ddptin=dptin+n(1)*n(1)
  need=ddptin+n(1)*n(1)
!
! Allocate some scratch space for this subroutine as one large array to
! be partitioned out.
!
  allocate(scr(need))
! 
! This array(s) is allocated locally and destroyed at the conclusion of the
! subroutine.  It is not available to the user.  If you wish to tinker with it, 
! do so at your own peril.
! 
! Calculate the functions and their derivatives using the scratch arrays.

  CALL twpt(scr(tin),scr(twtin),scr(ptin),scr(dptin),scr(ddptin),  &
            typwt,edge,n(1),n(1),prn(1))
  CALL nrmlze(scr(ptin),scr(dptin),scr(twtin),n(1))
  CALL tplmat(scr(ptin),scr(dptin),scr(twtin),scr(hin),n(1))
! 
! Move the scratch arrays to the saved space for the physical arrays.
!
  CALL tphy(scr(tin),scr(twtin),scr(ptin),scr(dptin),scr(ddptin),  &
            scr(hin),pt,pt(2),wt,f,df,ddf,h,n(1),nphy,prn(2))
!
! Deallocate the scratch space.
!
  deallocate(scr)
RETURN
END SUBROUTINE tlobato
















