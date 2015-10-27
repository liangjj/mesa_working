!deck fourier_basis.f90
!***begin prologue     fourier_basis
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           kinetic energy, hamiltonian, fourier basis
!***author             schneider, barry (nsf)
!***source
!***purpose            generate matrix elements
!***                   for a fourier expansion.
!***
!***description
!***references
!***routines called
!***end prologue       fourier_basis
  SUBROUTINE fourier_basis(pt,ke,h,v,eigv_0,eigvec_0,eigv,eigvec,    &
                           coord,n)
  USE dvr_global,     ONLY   : mass, iout,                           &
                               pi, one, two, three, box, nodiag
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: n    
  INTEGER                                :: i, j, k, l, m, delim, pre, ifac
  REAL*8, DIMENSION(n)                   :: pt, v, eigv_0, eigv
  REAL*8, DIMENSION(n,n)                 :: ke, h, eigvec_0, eigvec
  REAL*8                                 :: kii, kij, scale
  REAL*8                                 :: cosfac, sinfac_2
  REAL*8, DIMENSION(:), ALLOCATABLE      :: work
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=*)                      :: coord
  j=(n-1)/2
  kii=-pi*pi*(n*n-one)/(three*box*box)
  kij=-two*pi*pi/(box*box)
  DO  i=1,n
      ke(i,i)=kii
  END DO
  k=-j
  DO i=1,n
     l=-j
     DO m=1,i-1
        delim=k-l
        cosfac=cos(pi*delim/n) 
        sinfac_2=sin(pi*delim/n)
        sinfac_2=sinfac_2 * sinfac_2
        pre=1
        ifac=delim - 2*(delim/2)
        if(ifac == 1) THEN
           pre=-1
        END IF
        ke(i,m) = kij * pre * cosfac / sinfac_2
        l = l + 1
     END DO
     k = k + 1
  END DO
  DO i=1,n
     DO j=1,i
        ke(j,i) = ke(i,j)
     END DO
  END DO
  call vrmat(v,pt,coord,scale,n,1)
  ke=scale*ke
  IF(prn(3)) THEN
     title='normalized kinetic energy matrix for fourier basis' 
     CALL prntrm(title,ke,n,n,n,n,iout)
  END IF
  h=ke
  DO i=1,n
     h(i,i) = ke(i,i) + v(i)
  END DO
  IF(prn(3)) THEN
     title='normalized hamiltonian matrix for fourier basis' 
     CALL prntrm(title,h,n,n,n,n,iout)
  END IF
  IF(.not.nodiag) THEN
     ALLOCATE(work(5*n))
     call diarep(ke,eigv_0,eigvec_0,work,n)
     call diarep(h,eigv,eigvec,work,n)
     DEALLOCATE(work)
  END IF
END SUBROUTINE fourier_basis

