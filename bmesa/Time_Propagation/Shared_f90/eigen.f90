! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Eigen}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck eigen.f
!**begin prologue     eigen
!**date written       951229   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            eigenvalues and eigenvectors
!**
!**references
!**routines called
!**end prologue       eigen
  SUBROUTINE eigen(dim)
  USE io
  USE arnoldi_global
  IMPLICIT NONE
  INTEGER                                :: dim
  REAL*8, DIMENSION(:,:), ALLOCATABLE    :: htemp
  REAL*8, DIMENSION(:), ALLOCATABLE      :: scr
  INTEGER                                :: info
  ALLOCATE(scr(5*nphy(dim)))
!     diagonalize saving unitary transformation matrix in eigvec
  IF(row(dim) == 2) THEN
     ALLOCATE(htemp(row(dim),nphy(dim)))
     CALL cpy_3(grid(dim)%h,grid(dim)%e,htemp,nphy(dim))
     CALL dstev('v',nphy(dim),grid(dim)%e,htemp, &
                    grid(dim)%ev,nphy(dim),scr,info)
     DEALLOCATE(scr,htemp)
  ELSE IF(row(dim) == 3) THEN
     ALLOCATE(htemp(row(dim),nphy(dim)))
     htemp = grid(dim)%h
     CALL dsbev('v','l',nphy(dim),2,htemp,3,grid(dim)%e, &
                 grid(dim)%ev,nphy(dim),scr,info)
     DEALLOCATE(scr,htemp)
  ELSE IF(row(dim) == 4) THEN
     ALLOCATE(htemp(row(dim),nphy(dim)))
     htemp = grid(dim)%h
     CALL dsbev('v','l',nphy(dim),3,htemp,4,grid(dim)%e, &
                 grid(dim)%ev,nphy(dim),scr,info)
     DEALLOCATE(scr,htemp)
  ELSE
     grid(dim)%ev=grid(dim)%h
     title='hamiltonian in eigen'
     CALL prntrm(title,grid(dim)%ev,nphy(dim),nphy(dim), &
                                    nphy(dim),nphy(dim),iout)
     CALL dsyev('v','l',nphy(dim),grid(dim)%ev,nphy(dim), &
                                  grid(dim)%e,scr, &
                                  5*nphy(dim),info)
     DEALLOCATE(scr)
  END IF
  IF(info == 0) THEN
     IF(log_main(1)) THEN
        title='eigenvalues'
        CALL prntrm(title,grid(dim)%e,nphy(dim),1, &
                    nphy(dim),1,iout)
     END IF
     IF(log_main(2)) THEN
        title='eigenvectors'
        CALL prntrm(title,grid(dim)%ev,nphy(dim),nphy(dim), &
                    nphy(dim),nphy(dim),iout)
     END IF
  ELSE
     CALL lnkerr('error in diagonalization')
  END IF
END SUBROUTINE eigen




