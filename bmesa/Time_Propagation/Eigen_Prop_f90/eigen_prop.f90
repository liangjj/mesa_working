!documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{ Exact Propagator Code}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck eigen_prop
!**begin prologue     eigen_prop
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, time-propagation
!**
!**author             schneider, b. i.(nsf)
!**source             
!**purpose            time dependent schroedinger equation
!**                   eigenstate propagation with finite difference or
!**                   dvr space representation.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       eigen_prop
  SUBROUTINE eigen_prop
  USE arnoldi_global
  USE prop_global
  IMPLICIT NONE
  CHARACTER (LEN=16)                      :: fptoc
  CHARACTER (LEN=2)                       :: itoc
  INTEGER                                 :: i, bigs, bigv, tim
  INTEGER                                 :: iostat, length, len
  COMPLEX*16                              :: conjg
  REAL*8, DIMENSION(:), ALLOCATABLE       :: rtemp
!
!
  n3d=nphy(1)
  ALLOCATE(psi0(n3d),chi(n3d),soln_0(n3d),eig(n3d), &
!           u(n3d,n3d),work(n3d),tim_pts(ntreg+1), &
           work(n3d),tim_pts(ntreg+1), &
           auto_corr(ntreg),rtemp(max(n3d,ntreg+1)))
  stop
!  eig = grid(1)%eigv
!  u   = grid(1)%eigvec
!     Set up a file to hold the wavefunction for plotting
!
  call iosys('open plot as scratch',0,0,0,' ')
  call iosys('create real wavefunction on plot',ntreg*2*n3d,0,0,' ')
  DO tim=1,ntreg+1
     tim_pts(tim)=edge(tim)
  END DO
  IF(plot) then
     DO i=1,spdim
        write(iplot(1),*) grid(i)%pt
     END DO
     write(iplot(2),*) (tim_pts(i),i=1,ntreg+1)
  END IF
!
!     Begin the propagation.

  DO  tim=1,ntreg
      keywrd=itoc(tim)
      LEN=length(keywrd)
      keywrd='t'//keywrd(1:LEN)
      LEN=length(keywrd)
      keywrd='$v0('//keywrd(1:LEN)//')'
      t0=edge(tim)
      t1=edge(tim+1)
      deltat=t1-t0
      time(1)=secnds(0.0)
  
!  
!     Initialize the wavefunction at $t_{0}$ or 
!     read in its value from the disk at $t_{i-1}$.
!
      CALL cp_psi0(tim)
      time(2)=secnds(0.0)
      delta(1)=time(2)-time(1)
      CALL exp_prop('subtract',n3d)
      time(3)=secnds(0.0)  
      delta(2)=time(3)-time(2)
!  
!        compare the approximate and exact solutions where possible
  
      IF(i0stat == 'gaussian-pulse') THEN
         CALL moment
      END IF
      time(4)=secnds(0.0)  
      delta(3)=time(4)-time(3)
! Form the total solution and compute the autocorrelation
! function.
  
      CALL soln(tim,rtemp)
      time(5)=secnds(0.0)
      delta(4)=time(5)-time(4)
      call chk_nrm(chi,n3d)
      time(6)=secnds(0.0)
      delta(5)=time(6)-time(5)
      WRITE(iout,2) tim, (delta(i), i=1,5)

  END DO
!
!     End the propagation
!
  do i=1,ntreg
     rtemp(i) = auto_corr(i) * conjg(auto_corr(i))
!     rtemp(i) = - log(rtemp(i))/tim_pts(i+1)
  END DO
  write(iout,6) (rtemp(i),i=1,ntreg)
  CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
  IF(plot) then
     write(iplot(6),*) (rtemp(i),i=1,ntreg)
  END IF
!     Release the memory
!
  DEALLOCATE(psi0,chi,soln_0,eig, &
             u,work,tim_pts, &
             auto_corr,rtemp)

1 FORMAT('***********************************************'   &
         '*************************'                         &
         /,10X,'Time to pack the spatial Hamiltonian = '     &
               ,f15.8,/                                      &
         '***********************************************'   &
         '*************************')
2 FORMAT('***********************************************'     &
         '*************************'                           &
         /,10X,'Time Summary for Interval         = ',i4,      &
         /,10X,'time for right hand side          = ',f15.8,   &
         /,10X,'time to propagate                 = ',f15.8,/, &
         /,10X,'time to test wavefunction         = ',f15.8,/, &
         /,10X,'time to tabulate solution         = ',f15.8,/, &
         /,10X,'time to check normalization       = ',f15.8,/, &
         '***********************************************'     &
         '*************************')
3 FORMAT(5x,'points dimension = 'i2)
4 FORMAT(6e15.8)
5 FORMAT(a16)
6 FORMAT(/15x,'auto-correlation function = ',(/,5e15.8))
7 FORMAT(/,15x,'gamma = ',(/,5e15.8))
END SUBROUTINE eigen_prop
