!
!*deck fd_basis
!
!***begin prologue     fd_basis
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           finite difference
!***author             schneider, b. i.(nsf)
!***source             fdlib
!***purpose            calculate finite difference quantities
!***description        a three, five or seven point finite
!***                   difference is used to compute  
!                      the arrays are:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!            Array       Size                   Description
!            ____        ____                   ___________
!
!***          pt        nphy+1                  dvr points
!
!***          wt        nphy                    dvr weights
!
!***          ke       (row,nphy)               kinetic energy matrix
!
!***          eigv_0    nphy                    eigenvalues of kinetic energy matrix
!
!***          eigvec_0 (nphy,nphy)              eigenvectors of kinetic energy matrix
!
!***          h        (row,nphy)               hamiltonian matrix
!
!***          eigv      nphy                    eigenvalues of hamiltonian matrix
!
!***          eigvec_0 (nphy,nphy)              eigenvectors of hamiltonian matrix
!
!***          v         nphy                    diagonal matrix of potential
!
!***references

!***routines called    
!***end prologue       fd_basis
  Subroutine fd_basis(q_0,q,wt,ke,eigv_0,eigvec_0,h,eigv,eigvec,v, &
                      nphy,nglobal,row,coord)
  USE fd_global
  USE fd_prnt
  IMPLICIT NONE
  INTEGER                               :: nphy, nglobal, row
  CHARACTER(len=*)                      :: coord
  CHARACTER(len=80)                     :: title
  LOGICAL                               :: dollar, logkey
  REAL*8, DIMENSION(nphy)               :: v 
  REAL*8                                :: q_0
  REAL*8, DIMENSION(nglobal)            :: q, wt  
  REAL*8, DIMENSION(nphy)               :: eigv_0, eigv 
  REAL*8, DIMENSION(row,nphy)           :: ke, h
  REAL*8, DIMENSION(nphy,nphy)          :: eigvec_0, eigvec 
  REAL*8, DIMENSION(:), ALLOCATABLE     :: scr
  REAL*8, DIMENSION(:,:), ALLOCATABLE   :: htemp
  INTEGER                               :: info
  d=0.d0
  CALL points(q_0,q,wt,nphy)
  CALL vrmat(v,q,coord,dscale,nphy,1)
  IF(ndiff == 3) THEN
     CALL band3(ke,nphy)
  ELSE IF(ndiff == 5) THEN
     CALL band5(ke,nphy)
  ELSE IF(ndiff == 7) THEN
     CALL band7(ke,nphy)
  ELSE
     write(iout,1)
     stop
  END IF
  ke = dscale * ke
  h = ke
  h(1,:) = h(1,:) + v(:) 
  if(prn_fd_log(4)) then
     title='ke'
     call prntrm(title,ke,row,nphy,row,nphy,iout)
     title='h'
     call prntrm(title,h,row,nphy,row,nphy,iout)
  endif
  IF(.not.nodiag) then
     ALLOCATE(htemp(row,nphy),scr(5*nphy)) 
     IF(ndiff == 3) THEN
        CALL cpy_3(h,eigv,htemp,nphy)
        CALL dstev('v',nphy,eigv,htemp,eigvec,nphy,scr,info)
        CALL cpy_3(ke,eigv_0,htemp,nphy)
        CALL dstev('v',nphy,eigv_0,htemp,eigvec_0,nphy,scr,info)
     ELSE IF(ndiff == 5) THEN
        htemp = h
        CALL dsbev('v','l',nphy,2,htemp,3,eigv,eigvec,nphy,scr,info)
        htemp = ke
        CALL dsbev('v','l',nphy,2,htemp,3,eigv_0,eigvec_0,nphy,scr,info)
     ELSE IF(ndiff == 7) THEN
        htemp = h
        CALL dsbev('v','l',nphy,3,htemp,4,eigv,eigvec,nphy,scr,info)
        htemp = ke
        CALL dsbev('v','l',nphy,3,htemp,4,eigv_0,eigvec_0,nphy,scr,info)
     ELSE
        CALL lnkerr('error')
     ENDIF
     if(prn_fd_log(5)) then
        title='eigenvalues of h'
        call prntrm(title,eigv,nphy,1,nphy,1,iout)
        title='eigenvalues of ke'
        call prntrm(title,eigv_0,nphy,1,nphy,1,iout)
     endif
     if(prn_fd_log(6)) then
        title='eigenvectors of h'
        call prntrm(title,eigvec,nphy,nphy,nphy,nphy,iout)
        title='eigenvectors of ke'
        call prntrm(title,eigvec_0,nphy,nphy,nphy,nphy,iout)
     end if
     DEALLOCATE(htemp,scr)
  end if
1 format(/,1x,'Finite Difference Formula not Available')
END SUBROUTINE fd_basis

