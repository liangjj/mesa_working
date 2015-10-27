!
  MODULE fd_prnt
!***begin prologue     fd_prnt
!***date written       030603   (yymmdd)
!***revision date               (yymmdd)
!***keywords           finite difference
!***author             schneider, b. i.(nsf)
!***source             fdlib
!***purpose            global print variables for finite difference 
!***                   library
!***references
!***routines called    
!***end prologue       fd_prnt

  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!                   Keywords for Print Options
!
   CHARACTER (LEN=80), DIMENSION(7) :: prn_fd
   DATA prn_fd / 'global-points','global-weights','potential',  &
                 'hamiltonian','eigenvalues','eigenvectors','none' /
!
!
  CHARACTER (LEN=80), DIMENSION(7)      :: prn_fd_loc
  LOGICAL, DIMENSION(7)                 :: prn_fd_log
!
!
END MODULE fd_prnt
