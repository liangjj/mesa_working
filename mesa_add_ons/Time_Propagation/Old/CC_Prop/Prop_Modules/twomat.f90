lanczos_module.mod
!deck twomat.f
!***begin prologue     twomat
!***date written       030131   (yymmdd)
!***revision date               (yymmdd)
!***keywords           two electron radial integrals
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate two electron radial integrals
!***                   in a dvr representation by solving the
!***                   poisson equation for the inner integral
!***                   and then applying quadrature. 

!***references

!***routines called    iosys, util and mdutil
!***end prologue       twoel

  SUBROUTINE twomat(v_two,ke,ang_pot,pt,wt,nphy)
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: nphy
  REAL*8,   DIMENSION(nphy)              :: pt, wt, ang_pot
  REAL*8,   DIMENSION(nphy,nphy)         :: ke
  REAL*8,   DIMENSION(nphy_tri,0:l_max)  :: v_two
  REAL*8,   DIMENSION(:,:), ALLOCATABLE  :: t
  REAL*8,   DIMENSION(:,:), ALLOCATABLE  :: work
  REAL*8                                 :: scale, ptlst 
  INTEGER,  DIMENSION(:), ALLOCATABLE    :: ipvt
  INTEGER                                :: l, lval, i, j, info, count
!
  scale=-2.d0*mass
  ALLOCATE(t(nphy-1,nphy-1),work(nphy,5),ipvt(nphy))
  do l=0,l_max
!
!    put the second derivative sub-matrix into t and then scale it
!
     t(:,:)=ke(1:nphy-1,1:nphy-1)
     t = scale*t
     lval=l*(l+1)
!
!    add in the angular part
!
     do i=1,nphy-1
        t(i,i) = t(i,i) - lval*ang_pot(i)
     end do
!
!    invert t
!
     CALL dsytrf('l',nphy-1,t,nphy-1,ipvt,work,5*nphy,info)
     CALL dsytri('l',nphy-1,t,nphy-1,ipvt,work,info)
!
!    calculate the integrals
!
     work(:,1)=pt**l
     work(:,2)=1.d0/(pt*sqrt(wt))
     lval=2*l+1
     ptlst=1.d0/(pt(nphy)**lval)
     count=0
     do i=1,nphy-1
        do j=1,i
           count=count+1
           v_two(count,l) = -lval*t(i,j)*work(i,2)*work(j,2) + &
                             work(i,1)*work(j,1)*ptlst
        end do
     end do
     do j=1,nphy
        count=count+1
        v_two(count,l) = work(nphy,1)*work(j,1)*ptlst
     end do
  end do
  DEALLOCATE(t,work,ipvt)
END SUBROUTINE twomat
