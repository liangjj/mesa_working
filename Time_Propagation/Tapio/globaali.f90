!===========================================================================
! This module contains global variables including all parameters needed to 
! to run the code. The specific form of the potential must be changed in 
! subroutine "ve_prop" included in the file pot_propagation.f90. 
!===========================================================================

module globaali
  implicit none


  ! ===========  global parameters ==========================================================================================================
  character(LEN=5)                      :: indir      = 'DATA/'          ! input data directory
  character(LEN=5)                      :: outdir     = 'DATA/'          ! output data directory
  character(LEN=6)                      :: infile     = 'wf_abc'         ! input wavefunction data file 
  character(LEN=6)                      :: outfile    = 'wf_abc'         ! output wavefunction data file
  character(LEN=3), parameter           :: method     = 'SO2'            ! time propagation method (SO2) or (SO4)
  character(LEN=2), parameter           :: instate    = 'CC'             ! Options: (GA), (TF), (CC), (RN), (VT), (IN)
  real*8,     parameter                 :: C          = 100d0            ! nonlinearity constant
  real*8,     parameter                 :: lambdax    = 1d0              ! trap anisotropy omega_x / omega_0
  real*8,     parameter                 :: lambday    = 1d0              ! trap anisotropy omega_y / omega_0
  real*8,     parameter                 :: lambdaz    = 1d0              ! trap anisotropy omega_z / omega_0
  real*8,     parameter                 :: mutf       = 10               ! TF chemical potential (used only if instate='TF/VT')
  real*8,     parameter                 :: gamma      = 0d0              ! damping strength
  real*8,     parameter                 :: omega      = 0d0            ! angular velocity of the rotating frame
  real*8,     parameter                 :: tmax       = 4.5              ! end time of real time simulation (2*pi = 1 trap period)
  real*8,     parameter                 :: tol        = 0.007                ! tolerance for adaptive time stepping in imaginary time propagation
  real*8,     parameter                 :: maxiter    = 1e6              ! maximum number of iterations
  complex*16                            :: dt         = (0d0, -0.01d0)  ! time step: either pure real or pure imaginary
  complex*16, parameter                 :: iu         = (0d0, 1d0)       ! imaginary unit
  logical                               :: io         = .TRUE.           ! output on/off
  
  ! ================ MPI ====================================================================================================================
  integer                               :: rank, size, ierr              ! standard communicator variables   
  integer                               :: cart_comm, top, bot, lft, rgt ! variables for cartesian communicator
  integer                               :: root     = 0                  ! root process
  integer                               :: ndims    = 2                  ! number of dimensions
  integer, dimension(2)                 :: dims     = (/2,2/)            ! number of processes in each dimension
  logical                               :: reorder  = .true.             ! allow reordering of rows and columns 
  logical, dimension(2)                 :: periodic = (/0,0/)            ! say no to periodic boundary conditions 
  integer, dimension(128)               :: x_begind, x_endind, x_begblk, x_endblk    ! x-indexing variables: dimension >= dims(2)
  integer, dimension(128)               :: y_begind, y_endind, y_begblk, y_endblk    ! y-indexing variables: dimension >= dims(1)
  integer, dimension(2)                 :: coords                        ! process coordinates
  integer                               :: my_xind, my_yind              ! indices of "self process" 
  integer                               :: my_dimx, my_dimy              ! dimensions of "self process" 

  integer                               :: impdo                         ! dummy loop variable used in setting up the grid below 
  ! ========================  setup X grid ==================================================================================================
  integer,                   parameter  :: nregx       = 16              ! number of DVR blocks in x direction   
  integer,                   parameter  :: ppbx        = 7               ! number of points in each x-block
  real*8,                    parameter  :: Xbeg        = -8              ! -x boundary of the computational region
  real*8,                    parameter  :: Xend        =  8              ! +x boundary of the computational region
  integer, dimension(nregx), parameter  :: nptsx       = (/(ppbx,impdo=1,nregx)/) 
  real*8,  dimension(nregx+1)           :: boundariesx = (/( Xbeg + (impdo-1) * (Xend-Xbeg) / real(nregx) ,impdo=1,nregx+1)/)
  ! ========================  setup Y grid ==================================================================================================
  integer,                   parameter  :: nregy       = 16              ! number of DVR blocks in y direction
  integer,                   parameter  :: ppby        = 7               ! number of points in each y-block
  real*8,                    parameter  :: Ybeg        = -8             ! -y boundary of the computational region
  real*8,                    parameter  :: Yend        = 8              ! +y boundary of the computational region
  integer, dimension(nregy), parameter  :: nptsy       = (/(ppby,impdo=1,nregy)/) 
  real*8,  dimension(nregy+1)           :: boundariesy = (/( Ybeg + (impdo-1) * (Yend-Ybeg) / real(nregy) ,impdo=1,nregy+1)/)
  ! ========================  setup Z grid ==================================================================================================
  integer,                   parameter  :: nregz       = 16             ! number of DVR blocks in z direction
  integer,                   parameter  :: ppbz        = 7               ! number of points in each z-block
  real*8,                    parameter  :: Zbeg        = -4             ! -z boundary of the computational region
  real*8,                    parameter  :: Zend        =  4             ! +z boundary of the computational region
  integer, dimension(nregz), parameter  :: nptsz       = (/(ppbz,impdo=1,nregz)/)
  real*8,  dimension(nregz+1)           :: boundariesz = (/( Zbeg + (impdo-1) * (Zend-Zbeg) / real(nregz) ,impdo=1,nregz+1)/)


  ! ============ global variables ===========================================================================================================
  real*8,  dimension(2), parameter      :: bridgepts  = (/-1d0,1d0/)     ! bridge (end) points
  real*8,  dimension(3)                 :: pq                            ! parameters needed for 4th order time propagation
  real*8                                :: myy, err, errold              ! chemical potential and its error
  complex*16                            :: aika, ang                     ! aika = current time in simulation, ang = <angular momentum> 
  integer                               :: dimx, dimy, dimz              ! extent of each spatial dimension

  ! =========== global vectors and matrices =================================================================================================
  real*8,    allocatable, dimension(:,:)      :: pwx                     ! points and weights (pt,wt) bridges removed
  real*8,    allocatable, dimension(:,:)      :: pwy                     ! points and weights (pt,wt) bridges removed
  real*8,    allocatable, dimension(:,:)      :: pwz                     ! points and weights (pt,wt) bridges removed
  complex*16,allocatable, dimension(:,:,:)    :: wf, hwf, hpsi,locerr, ehwf     ! wf = wavefunction ( the rest are derived from it)
  complex*16,allocatable, dimension(:,:,:)    :: vext                    ! external potential such as the trap 
  complex*16,allocatable, dimension(:,:,:)    :: vprop                   ! external potential propagator
  
 
  ! =========== derived type for block data ================================================================================================
  type block
     real*8,    allocatable, dimension(:)     :: dummy, pt, wt          ! block points and weights
     real*8,    allocatable, dimension(:,:)   :: p, dp, ddp             ! DVR functions and derivatives 
     real*8,    allocatable, dimension(:,:)   :: ke                     ! kinetic energy block matrices 
     real*8,    allocatable, dimension(:,:)   :: dx, dy, dz             ! momentum operator block matrices    
     complex*16,allocatable, dimension(:,:)   :: mprop                  ! momentum block propagators
     complex*16,allocatable, dimension(:,:)   :: mpropq                 ! momentum block propagators
     complex*16,allocatable, dimension(:,:)   :: kprop                  ! kinetic energy block propagators
     complex*16,allocatable, dimension(:,:)   :: kpropq                 ! kinetic energy block propagators
  end type block
  ! ========== derived type for angular momentum block data ================================================================================
  type blocks
     type(block), allocatable, dimension(:)   :: blk
  end type blocks
  ! =========== x,y, and z blocks ==========================================================================================================
  type(block),  allocatable, dimension(:)     :: lohkox
  type(block),  allocatable, dimension(:)     :: lohkoy
  type(block),  allocatable, dimension(:)     :: lohkoz
  ! ============ angular momentum blocks ===================================================================================================  
  type(blocks), allocatable, dimension(:)     :: lxprop
  type(blocks), allocatable, dimension(:)     :: lyprop

end module globaali

