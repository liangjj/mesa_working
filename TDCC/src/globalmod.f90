module globalmod
  USE nrtype
  implicit NONE

  !----constants---------------------------- 

  INTEGER(I4B), parameter :: cndim = 2 ! number of dimensions for MPI cartesian communicator
  INTEGER(I4B), parameter :: ndim = 2  ! number of spatial physical dimensions (i.e. without angular part)

  REAL(DP) , parameter :: zero=0.d0, one=1.d0, &
       two=2.d0, half=0.5d0, quarter=0.25d0, &
       convfs=0.02418884d0, conve=27.211608d0, &
       convl=0.529177249d-10

  !------------for cartesian communicator comm2d-------------

  INTEGER(I4B) :: coords(cndim), cdims(cndim), nbrup, nbrdown, &
       nbrleft, nbrright, ncpux, ncpuy, comm2d

  INTEGER(I4B) :: cperiods(cndim), creorder 

  !------------propagation parameters-------------

  REAL(DP) ::   intensity_x, omegax, durationx, phasex, &
       intensity_y, omegay, durationy, phasey, &
       Z1, alpha, RR, VRi, alphai,tstart, tstop, &
       wavelx, wavely, efieldx, efieldy, TCx, TCy, &
       dt, efldx, efldy, gfactor, xexp, yexp, zexp, &
       pote, kinex, kiney, kinez, potee

  INTEGER(I4B) :: nabx, naby, nabz, ntim, ntrec, ntrep, ntprob,  &
       ncount, nrcount, ntemp1, ntemp2, nelem, &
       nelemx, nelemy, nelemz, it, nelemyz, nelemxy, npointxy, &
       npointx, npointy, npointz, npointyz, npointxz, ix, iy                

  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE :: bdstart, bdstop, bdslope
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE, SAVE :: numfunc, ntarget, nst1, nst2


  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE :: prob1d, probx1d
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, SAVE ::  prob2d, probx2d, proby2d
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE ::  prob3d, probx3d, proby3d, probz3d

  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE :: prob1dtemp, probdens1d
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, SAVE :: prob2dtemp, probdens2d 
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, SAVE :: prob2dtempxy, probdens2dxy 


  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE :: cprob1d 
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, SAVE :: cprob2d, cprobpx2d, cprobpy2d
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, SAVE ::  psitemp2d, psitemp2dxy, psitemp2dyz, psitemp2dxz
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: cprob3d, cprobpx3d, cprobpy3d, cprobpz3d    


  REAL(DP) ::  ctemp_mat0(200,200), ctemp_mat1(200,200), &
       ctemp_mat2(200,200), ctemp_mat3(5000), &
       ctemp_mat4(200,200,200),ctemp_mat5(200,200,200),ctemp_mat6(200,200,200)

  REAL(DP)  ::   FV1(500), FV2(500)

  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE ::  psitemp1d

  !----------Boundary condition at r=0 or not------------------------------

  logical :: r0BC, nonlinear, BEC, Strong_Field


  !-----------Complete spatial grid and mapping (for "ndim" dimensions)------------------

  REAL(DP) , dimension(:), ALLOCATABLE :: ax, ay, dipx, dipy, bmat_eigval
  REAL(DP) , dimension(:,:), ALLOCATABLE :: grid , bmat_eigvec, bmat
  REAL(DP) , dimension(:,:), ALLOCATABLE :: factor 
  INTEGER(I4B), dimension(:,:), ALLOCATABLE :: point_order
  INTEGER(I4B), dimension(:), ALLOCATABLE :: regstart, regstop, ngridstart, ngridstop, ngstart, ngstop


  !-------------Total number of grid points for each dimension--------------- 

  integer , dimension(:), ALLOCATABLE :: ntot


  !----------Kinetic-energy matrix, eigenvector, eigenvalue for--------------- 
  !----------each dinmension, each region: (ndim,ith-region,nfun,nfun)---------------


  REAL(DP), dimension(:,:,:,:), ALLOCATABLE :: Kin_eigvec, Kin_eigvec_inv

  REAL(DP), dimension(:,:,:,:), ALLOCATABLE :: Kin_operator_dt, &
       Kin_operator_halfdt, &
       Kin_operator_quarterdt

  REAL(DP), dimension(:,:,:), ALLOCATABLE :: Kin_eigval_dt, &
       Kin_eigval_halfdt, &
       Kin_eigval_quarterdt

  !--------------wavefunction-------------------------

  REAL(DP), dimension(:,:,:), ALLOCATABLE  :: psi3d 
  REAL(DP), dimension(:,:,:), ALLOCATABLE  :: psi03d 
  REAL(DP), dimension(:,:,:), ALLOCATABLE  :: coeff3d
  REAL(DP), dimension(:,:,:), ALLOCATABLE  :: wholepsi3d


  !--------------potential-------------------------

  REAL(DP), dimension(:,:,:), ALLOCATABLE  ::   cpot3d

  REAL(DP), dimension(:,:,:), ALLOCATABLE  ::   rpot3d

  !--------------mapping index and arrays-------------------------


  INTEGER(I4B), dimension(:,:,:), ALLOCATABLE :: nindex



  !--------------Parameters for constructing grid-------------------------

  INTEGER(I4B)  :: start,start2,end


  !--------------for two electron interactions-------------------------


  INTEGER(I4B), DIMENSION(:), ALLOCATABLE, SAVE :: l1mat,l2mat,ltotmat , ivv1   

  INTEGER(I4B), SAVE :: LL,l0,l1,l2, Lp,l1p, l2p, nwave, nnp, mm, nparity, &
       ir1, ir2, starti, startj

  REAL(DP), SAVE  ::   rmin, rmax

  REAL(DP), dimension(:), ALLOCATABLE  ::  eiger, eigei, fvv1, vv1i, vv2i

  REAL(DP), dimension(:,:), ALLOCATABLE  ::   amat, BB, BT, F1, F2

  REAL(DP), dimension(:,:,:), ALLOCATABLE  ::  EH, F111, F222, F0

  REAL(DP), dimension(:,:,:,:), ALLOCATABLE  ::   VL, GG, GT, HH, HT

  REAL(DP), dimension(:,:,:,:), ALLOCATABLE  ::  CGG

  REAL(DP), dimension(:,:,:), ALLOCATABLE  ::   EG

  REAL(DP), dimension(:,:), ALLOCATABLE  ::   camat

  !----------------------------------------------------------------------


end module globalmod
