!
!*deck atomic_basis
!
!***begin prologue     atomic_basis
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           time, finite element dvr, orthogonal polynomial
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            calculate dvr functions, derivatives and matrix elements
!***description        this routine computes the dvr arrays
!***                   needed to discretize a one-dimensional
!***                   hamiltonian using the finite element dvr
!***                   method.  the number of unknowns, nphy,
!***                   is not the unique number of points in the grid
!***                   but what is left after accounting for boundary
!***                   conditions which constrain the wavefunction. 
!                      the arrays are:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!            Array       Size                   Description
!            ____        ____                   ___________
!***          pt_0        1                     leftmost dvr point
!
!***          pt_n        1                     rightmost dvr point
!
!***          pt        nphy                    dvr points
!
!***          wt        nphy                    dvr weights
!
!***          f        (nphy,nphy)              dvr functions
!
!***          df       (nphy,nphy)              first derivative of dvr functions
!
!***          ddf      (nphy,nphy)              second derivative of dvr functions
!
!***          ke       (nphy,nphy)              FEDVR kinetic energy matrix. 
!                                               Contribution from Bloch
!                                               operator added.
!
!***          p_mom    (nphy,nphy)              FEDVR First derivative matrix.
!                                               Contribution with Bloch 
!                                               operator added.  Needed for 
!                                               interactions with external 
!                                               EM field.  
!
!***          cou_pot   nphy                    diagonal matrix of 1/r
!
!***          ang_pot   nphy                    diagonal matrix of 1/(r*r)
!
!***          srf_prm   2                       primitive surface values of dvr basis
!
!***          hmat(nphy,nphy,0:l_orb_max)       one electron integrals
!
!***          v_twoel(nphy*(nphy+1)/2,0:l_max)  radial two-electron integrals
!***references

!***routines called    
!***end prologue       atomic_basis
  Subroutine atomic_basis(nphy,nglobal)
  USE dvr_global
  IMPLICIT NONE
!
!                   The Main Arrays Returned to the User
!
  INTEGER                               :: nphy, nglobal
  REAL*8                                :: pt_0, pt_n
  REAL*8, DIMENSION(:),     ALLOCATABLE :: pt
  REAL*8, DIMENSION(:),     ALLOCATABLE :: wt, cou_pot, ang_pot
  REAL*8, DIMENSION(:),     ALLOCATABLE :: rho, v_poisson
  REAL*8, DIMENSION(:,:),   ALLOCATABLE :: v_two
  REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: h_one
  REAL*8, DIMENSION(:,:),   ALLOCATABLE :: f, df, ddf, ke, p_mom
  REAL*8, DIMENSION(2)                  :: srf_prm
  INTEGER                               :: len, lenth
  ALLOCATE(pt(nphy),wt(nphy),cou_pot(nphy),ang_pot(nphy))
  ALLOCATE(f(nphy,nphy),df(nphy,nphy),ddf(nphy,nphy),ke(nphy,nphy), &
           p_mom(nphy,nphy))
  CALL Onel(pt_0,pt_n,pt,wt,f,df,ddf,ke,p_mom,cou_pot,ang_pot,srf_prm, &
            nphy,nglobal)
  ALLOCATE(h_one(nphy,nphy,0:l_orb_max))
  CALL Onemat(h_one,ke,p_mom,cou_pot,ang_pot,nphy)
  if(drctv == 'all-integrals') then
     ALLOCATE(v_two(nphy_tri,0:l_max))
     CALL Twomat(v_two,ke,ang_pot,pt(2),wt,nphy)
     DEALLOCATE(v_two)
  ELSE IF(drctv == 'poisson-equation') then
     CALL poisson(ke,ang_pot,pt(2),f,wt,nphy)
  end if
  DEALLOCATE(h_one,f,df,ddf,ke,p_mom,pt,wt,cou_pot,ang_pot)
END SUBROUTINE atomic_basis
