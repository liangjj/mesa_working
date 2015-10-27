!deck phy.f
!***begin prologue     phy
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            physical hamiltonian and potential.
!***
!***references

!***routines called
!***end prologue       phy
  SUBROUTINE phy(x,xwt,px,dpx,ddpx,k_mat,v_coul,v_ang,pt_0,pt, &
                 wt,f,df,ddf,ke,cou_pot,ang_pot,srf_prm,nglobal, &
                 nphy,start)
  IMPLICIT NONE
  INTEGER                             :: nglobal, nphy, ndel
  REAL*8, DIMENSION(nglobal)          :: x, xwt
  REAL*8, DIMENSION(nglobal,nglobal)  :: px, dpx, ddpx, k_mat
  REAL*8, DIMENSION(nglobal)          :: v_coul, v_ang
  REAL*8                              :: pt_0
  REAL*8, DIMENSION(nphy)             :: pt, wt
  REAL*8, DIMENSION(nphy,nphy)        :: f, df, ddf, ke
  REAL*8, DIMENSION(nphy)             :: cou_pot, ang_pot
  REAL*8, DIMENSION(2)                :: srf_prm
  INTEGER                             :: start
  CHARACTER (LEN=80)                  :: title
  pt_0=x(1)
  pt(1:nphy)=x(start:nglobal)
  wt(1:nphy)=xwt(start:nglobal)
  f(1:nphy,1:nphy)=px(start:nglobal,start:nglobal)
  df(1:nphy,1:nphy)=dpx(start:nglobal,start:nglobal)
  ddf(1:nphy,1:nphy)=ddpx(start:nglobal,start:nglobal)
  ke(1:nphy,1:nphy)=k_mat(start:nglobal,start:nglobal)
  cou_pot(1:nphy)=v_coul(start:nglobal)
  ang_pot(1:nphy)=v_ang(start:nglobal)
  srf_prm(1)=px(1,1)
  srf_prm(2)=px(nglobal,nglobal)
END SUBROUTINE phy



