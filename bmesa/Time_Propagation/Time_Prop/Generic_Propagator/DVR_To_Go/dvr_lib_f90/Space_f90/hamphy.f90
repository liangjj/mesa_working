!deck hamphy.f
!***begin prologue     hamphy
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            physical hamiltonian and potential.
!***
!***references

!***routines called
!***end prologue       hamphy
  SUBROUTINE hamphy(x,xwt,px,dpx,ddpx,hmat,vmat,kmat,pmat, &
                   pt_0,pt_n,pt,wt,f,df,ddf,p_mom,h,v,ke,srf_prm, &
                   nglobal,nphy,start)
  IMPLICIT NONE
  INTEGER                             :: nglobal, nphy, ndel
  REAL*8, DIMENSION(nglobal)          :: x, xwt
  REAL*8, DIMENSION(nglobal,nglobal)  :: px, dpx, ddpx, hmat, kmat, pmat
  REAL*8, DIMENSION(nglobal)          :: vmat
  REAL*8                              :: pt_0, pt_n
  REAL*8, DIMENSION(nphy)             :: pt, wt
  REAL*8, DIMENSION(nphy,nphy)        :: f, df, ddf, p_mom, h, ke
  REAL*8, DIMENSION(nphy)             :: v
  REAL*8, DIMENSION(2)                :: srf_prm
  INTEGER                             :: start
  CHARACTER (LEN=80)                  :: title
  pt_0=x(1)
  pt_n=x(nglobal)
!  CALL copy(x(start),pt,nphy)
  pt(1:nphy)=x(start:nglobal)
!  CALL copy(xwt(start),wt,nphy)
  wt(1:nphy)=xwt(start:nglobal)
!  CALL mmove(px(start,start),f,nphy,nphy,nglobal,nphy)
  f(1:nphy,1:nphy)=px(start:nglobal,start:nglobal)
!  CALL mmove(dpx(start,start),df,nphy,nphy,nglobal,nphy)
  df(1:nphy,1:nphy)=dpx(start:nglobal,start:nglobal)
!  CALL mmove(ddpx(start,start),ddf,nphy,nphy,nglobal,nphy)
  ddf(1:nphy,1:nphy)=ddpx(start:nglobal,start:nglobal)
!  CALL mmove(hmat(start,start),h,nphy,nphy,nglobal,nphy)
  h(1:nphy,1:nphy)=hmat(start:nglobal,start:nglobal)
  ke(1:nphy,1:nphy)=kmat(start:nglobal,start:nglobal)
  p_mom(1:nphy,1:nphy) = pmat(start:nglobal,start:nglobal)
!  CALL copy(vmat(start),v,nphy)
  v(1:nphy)=vmat(start:nglobal)
  srf_prm(1)=px(1,1)
  srf_prm(2)=px(nglobal,nglobal)
END SUBROUTINE hamphy



