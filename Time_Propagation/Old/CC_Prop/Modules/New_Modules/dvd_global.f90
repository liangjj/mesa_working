!
MODULE dvd_global
!***begin prologue     dvd_global
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           davidson, eigenvalues
!***author             schneider, b. i.(nsf)
!***source             dvdlib
!***purpose            global variables for davidson library
!***description        this routine defines the global variables
!***                   and data needed for the dvdlib
!
!***references

!***routines called    
!***end prologue       dvd_global
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!***     variable      type       size                 description
!***     _______       ____       ____             ___________
!***      h             real       nphy           one body hamiltonian in
!***                                              i coordinate with diagonal
!***                                              elements set to zero
!***      diag          real       n              diagonal part of matrix
!***      vec           real    n*maxvec          davidson vectors
!***      hvec          real    n*maxvec          action of the hamiltonian 
!                                                 on pvec
!***      b, bwrk       real    maxvec*maxvec     small hamiltonian matrix 
!                                                 and a copy
!***      eigwrk        real    maxvec            eigenvalues of small matrix
!***      work          real    max(5*maxvec,     scratch array
!***                                n_dvd*maxvec)
!***      svec          real    max(maxvec*maxvec small scratch matrix
!***                    real        n*maxvec)
!***      resid         real    n_dvd*maxvec      residual vectors
!***      eig           real    nroot             converged eigenvalues
!***      cnverg        real                      convergence criterion    
!                                                 for a root
!***      thresh        real                      overlap tolerance for     
!                                                 accepting a new davidson  
!                                                 vector
!***      n_dvd         integer                   matrix size
!***      nphy          integer ndim              dimension of h matrices
!***      nroot         integer                   number of roots to get
!***      ntrial        integer                   number of trial vectors
!***      nattim        integer                   number of roots to        
!                                                 converge at a pass
!***      maxit         integer                   maximum size of davidson  
!                                                 space
!***      maxvec        integer                   maximum number of vectors 
!***      nblck         integer                   preconditioner maximum    
!                                                 block size available for  
!                                                 storage
!***      prn_dvd       logical                   print flags   
!***
!***      precon        character                 type of preconditioner
!
  INTEGER, PARAMETER                         :: lenbuf = 1000000
  CHARACTER (LEN=1600)                       :: card
  CHARACTER (LEN=80)                         :: cpass, precon
  CHARACTER (LEN=80)                         :: title
  INTEGER                                    :: n_dvd, dim
  INTEGER                                    :: ntrip, nblck, nroot, &
                                                maxvec, maxit,       &
                                                nattim, ntrial
!
  INTEGER                                    :: npass, nleft, rtdone, root0, rootn
  INTEGER                                    :: num2do, nin, used, ipass, iter
  INTEGER                                    :: left, n2rd, numnew, addvec, chkno
  INTEGER                                    :: begin, size, nout, con, uncon, min2kp
  INTEGER                                    :: errno
  CHARACTER*8, DIMENSION(3)                  :: code
  CHARACTER*16                               :: drctv
  LOGICAL                                    :: orth
  REAL*8                                     :: rep
!
  REAL*8                                     :: cnverg, thresh
!
  REAL*8,  DIMENSION(:),       ALLOCATABLE   :: diag, eig, eigwrk
  REAL*8,  DIMENSION(:,:),     ALLOCATABLE   :: b, bwrk, work, svec
  REAL*8,  DIMENSION(:,:),     ALLOCATABLE   :: vec, hvec, resid
  REAL*8, DIMENSION(:),        ALLOCATABLE   :: scratch  
  INTEGER, DIMENSION(:,:),     ALLOCATABLE   :: ind  
  REAL*8                                     :: error, maxerr
!
!
  data code / 'trial:', 'dvr-ene:', 'dvr-vec:' / 
END MODULE dvd_global
