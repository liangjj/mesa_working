! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{ R_mat Code}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck r_mat
!**begin prologue     r_mat
!**date written       030525   (yymmdd)
!**revision date               (yymmdd)
!**keywords           r-matrix, dvr, fd
!**
!**author             schneider, b. i.(nsf)
!**source             R_Matrix_f90
!**purpose            Calculate R-matrix and related information.
!**description
!**references
!**routines called    iosys, util and mdutil
!**end prologue       r_mat
  Subroutine r_mat(nume)
  USE r_matrix_global
  USE dvr_shared
  IMPLICIT NONE
  INTEGER                          :: nume
  INTEGER                          :: i, mul, tst
  rmat=0.d0
  DO i=1,nphy(1)
     rmat = rmat + grid(1)%srf(i,2) * grid(1)%srf(i,2)/( grid(1)%eigv(i) - energy(nume))
  END DO
  rmat = .5d0*rmat
  k=sqrt(2.d0*energy(nume))
  rho = k * grid(1)%pt(nphy(1))
!
! rbes returns the ricatti-bessel functions with the asymptotic
! form given in the NBS handbook.  Note that this requires a sign
! change for irregular function depending on the l value.
!
  call rc1bes(rho,jbes,djbes,ybes,dybes,angmom,ltop, &
              'derivatives',.false.)
!  write(iout,*) jbes
!  write(iout,*) djbes
!  write(iout,*) ybes
!  write(iout,*) dybes
!
! Fix the sign of the irregular function so it goes as cos(kr -l*pi/2)
! for all l.
!
  mul=-1
  tst = angmom - 2*(angmom/2)
  if (tst == 1) then
      mul=1
  endif
!
! Get the derivatives wrt r not rho.
!
  ybes(angmom)=mul*ybes(angmom)
  dybes(angmom)=mul*dybes(angmom)
  djbes(angmom)=k*djbes(angmom)
  dybes(angmom)=k*dybes(angmom)
  kmat(nume)=( rmat*djbes(angmom) - jbes(angmom) ) / &
             ( ybes(angmom) - rmat*dybes(angmom) )
  phase(nume)=atan(kmat(nume))
  write(iout,1) energy(nume), kmat(nume), phase(nume)
1 format(/,5x,'energy = ',e15.8,1x,'k-matrix = ',e15.8,1x,'phase shift= ',e15.8)  
END SUBROUTINE r_mat
