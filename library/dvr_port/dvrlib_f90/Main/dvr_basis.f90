!
!*deck dvr_basis
!
!***begin prologue     dvr_basis
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
!                      Array       Size         Description
!                      ____        ____         ___________
!
!***                    pt_0       1            leftmost dvr point
! 
!***                    pt        nphy          dvr points
!
!***                    wt        nphy          dvr weights
!
!***                    f        (nphy,nphy)    dvr functions
!
!***                    df       (nphy,nphy)    first derivative of dvr functions
!
!***                    ddf      (nphy,nphy)    second derivative of dvr functions
!
!***                    ke       (nphy,nphy)    FEDVR kinetic energy matrix.
!                                               Contribution from Bloch operator
!                                               added.
!
!***                    p_mom    (nphy,nphy)    FEDVR first derivative matrix.
!                                               Contribution of Bloch operator
!                                               added.  
!
!***                    eigv_0    nphy          eigenvalues of ke if requested
!
!***                    eigvec_0 (nphy,nphy)    eigenvectors of ke if requested
!
!***                    h        (nphy,nphy)    FEDVR h = ( ke + v )
!
!***                    eigv      nphy          eigenvalues of h if requested
!
!***                    eigvec   (nphy,nphy)    eigenvectors of h if requested
!
!***                    v         nphy          diagonal dvr potential matrix
!
!***                    srf_prm   2             primitive surface values of dvr basis
!
!***                    srf_0    (2,nphy)       surface values of ke eigenvectors
!
!***                    srf      (2,nphy)       surface values of h eigenvectors
!
!***                   If we are dealing with the time coordinate, some
!***                   of these arrays are not relevant.
!***references

!***routines called    
!***end prologue       dvr_basis
  Subroutine dvr_basis(pt_0,pt,wt,f,df,ddf,ke,p_mom,eigv_0,eigvec_0, &
                       h,eigv,eigvec,v,srf_prm,srf_0,srf, &
                       coord,nphy,nglobal)
  USE dvr_global
  IMPLICIT NONE
!
!                   The Main Arrays Returned to the User
!
  INTEGER                        :: nphy, nglobal
  REAL*8                         :: pt_0
  REAL*8, DIMENSION(nphy)        :: pt
  REAL*8, DIMENSION(nphy)        :: wt, eigv_0, eigv, v
  REAL*8, DIMENSION(nphy,nphy)   :: f, df, ddf, ke, p_mom, eigvec_0, h, eigvec
  REAL*8, DIMENSION(2)           :: srf_prm
  REAL*8, DIMENSION(2,nphy)      :: srf_0, srf
  INTEGER                        :: len, chrlen
  CHARACTER(LEN=*)               :: coord
  len=chrlen(coord)
  IF(coord(1:1) == 't') THEN
     CALL tlobato(pt,wt,f,df,ddf,h,nphy)
  ELSE
     CALL lobatto(pt_0,pt,wt,f,df,ddf,ke,p_mom,eigv_0,eigvec_0, &
                     h,eigv,eigvec,v,srf_prm,srf_0,srf, &
                     coord(1:len),nphy,nglobal)
  END IF
  RETURN
END SUBROUTINE dvr_basis
