!deck poisson.f
!***begin prologue     poisson
!***date written       030131   (yymmdd)
!***revision date               (yymmdd)
!***keywords           two electron radial integrals
!***author             schneider, b. i.(nsf)
!***source
!***purpose            solve poisson equation using dvr functions

!***references

!***routines called    iosys, util and mdutil
!***end prologue       twoel

  SUBROUTINE poisson(ke,ang_pot,f,pt,wt,nphy)
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: nphy
  REAL*8,   DIMENSION(nphy)              :: pt, wt, ang_pot
  REAL*8,   DIMENSION(nphy,nphy)         :: ke, f
  REAL*8,   DIMENSION(:), ALLOCATABLE    :: rho, v
  REAL*8,   DIMENSION(:,:), ALLOCATABLE  :: t
  REAL*8,   DIMENSION(:,:), ALLOCATABLE  :: work
  REAL*8                                 :: scale, ptlst, a 
  INTEGER,  DIMENSION(:), ALLOCATABLE    :: ipvt
  INTEGER                                :: l, l_barrier, i, info
  INTEGER                                :: l_fac
  CHARACTER(LEN=80)                      :: title
!
  ALLOCATE(rho(nphy),v(nphy))
  CALL rhofun(rho,pt,nphy)
  scale=-2.d0*mass
  ALLOCATE(t(nphy-1,nphy-1),work(nphy,5),ipvt(nphy))
  do l=0,l_max
!
!    put the second derivative sub-matrix into t and then scale it
!
     t(:,:)=ke(1:nphy-1,1:nphy-1)
     t = scale*t
     l_barrier=l*(l+1)
!
!    add in the angular part
!
     do i=1,nphy-1
        t(i,i) = t(i,i) - l_barrier*ang_pot(i)
     end do
!
!    solve linear equations
!
     CALL dsytrf('l',nphy-1,t,nphy-1,ipvt,work,5*nphy,info)
     CALL dsytrs('l',nphy-1,1,t,nphy-1,ipvt,rho,1,info)
!
!    calculate the potential
!
     l_fac=2*l+1
     work(:,1)=pt**(l+1)
     ptlst=1.d0/(pt(nphy)**l_fac)
     a=0.d0
     do i=1,nphy
        a = a + rho(i)*wt(i)*work(i,1)
     end do
     a=a/ptlst
     v(nphy)= a*work(nphy,1)
     do i=1,nphy-1
        v(i) = rho(i)*f(i,i) + a*work(i,1)
     end do
     v = v / pt 
     title='potential'
     call prntrm(title,v,nphy,1,nphy,1,iout)
  end do
  DEALLOCATE(t,work,ipvt)
  DEALLOCATE(rho,v)
END SUBROUTINE poisson
