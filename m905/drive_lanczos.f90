! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Main Program for Testing Lanczos Code}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck drive_lanczos
!**begin prologue     drive_lanczos
!**date written       030824   (yymmdd)
!**revision date               (yymmdd)
!**keywords           
!**author             schneider, b. i.(nsf)
!**source             lanczos
!**purpose            driver for lanczos diagonalization
!**description
!**references
!**routines called    iosys, util and mdutil
!**end prologue       drive_lanczos
  PROGRAM drive_lanczos
  USE lanczos_global
  IMPLICIT NONE
  CHARACTER (LEN=1)                         :: itoc
  LOGICAL                                   :: dollar
  INTEGER                                   :: intkey, i, j, ierr
  CHARACTER*4096                            :: ops 
  CALL Drum
  CALL IOsys('read character options from rwf',-1,0,0,ops)
  WRITE(iout,*)
  n=intkey(ops,'size-of-matrix',2,' ')
  iter=intkey(ops,'number-of-lanczos-iterations',n,' ')
  write(iout,1) n, iter    
!
  ALLOCATE(v(n,0:iter),a(iter),b(iter),eig(iter),vec(iter,iter),hvec(n), &
           scr(n),matrix(n,n))
  CALL rdmat
  CALL init_vec
  DO i=1,iter
     CALL lanczos(i)
     vec = 0.d0
     DO j=1,i
        vec(j,j) = 1.d0
     END DO
     CALL imtql2(iter,i,a,b,vec,ierr)
     title='eigenvalues of lanczos polynomial'
     CALL prntrm(title,a,i,1,i,1,iout)
     title='eigenvectors of lanczos polynomial'
     call prntrm(title,vec,i,i,iter,iter,iout)
  END DO
  call chainx(0)
  stop
1    FORMAT(/,20X,'Lanczos Code', &
            /,10x,'size of matrix               = ',i3, &
            /,10x,'maximum number of iterations = ',i3)
END PROGRAM drive_lanczos





