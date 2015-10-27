program GPfemdvr
  !=========================================
  ! Driver for solving 3D NLSE using FEM DVR
  ! Author: TPS 
  ! - parallel version, Jan. 2007 
  !=========================================
  use globaali, only: outdir, outfile, io
  use globaali, only: nptsx, nptsy, nptsz 
  use globaali, only: nregx, nregy, nregz 
  use globaali, only: wf, vext, vprop 
  use globaali, only: lohkox, lohkoy, lohkoz
  use globaali, only: lxprop, lyprop
  use globaali, only: pwx, pwy, pwz, dimx, dimy, dimz
  use globaali, only: C, omega, lambday, lambdaz, instate, lxprop, lyprop
  ! ===== MPI=====
  use globaali, only: root, rank, size, ierr                            
  use globaali, only: cart_comm, ndims, dims, periodic, reorder, coords 
  use globaali, only: top, bot, lft, rgt                                
  use globaali, only: x_begind, x_endind, y_begind, y_endind            
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk            
  use globaali, only: my_xind, my_yind 
  use globaali, only: my_dimx, my_dimy      
  implicit none

  real*8  :: tbeg, tend  ! for timing the CPU usage 
  integer :: dumind      ! dummy variable

  !======================================
  ! Initialize MPI 
  !======================================
  include 'mpif.h'                                              
  call MPI_Init(ierr)                                           
  call MPI_Comm_Rank(MPI_Comm_World, rank, ierr)               
  call MPI_Comm_Size(MPI_Comm_World, size, ierr)              
  call MPI_Cart_create(MPI_Comm_World, ndims, dims, periodic, reorder, cart_comm, ierr) 
  call MPI_Cart_Coords(cart_comm, rank, ndims, coords, ierr)
  call MPI_Cart_shift (cart_comm, 0, 1, top, bot, ierr )        
  call MPI_Cart_shift (cart_comm, 1, 1, lft, rgt, ierr )  
  my_xind = coords(2) + 1
  my_yind = coords(1) + 1
  print*, 'Rank',rank, 'has index (x,y):', my_xind, my_yind
!  print*, 'Process:', rank, 'is alive.'                         
!  print*, 'Process:', rank ,' has top:',top, 'and bottom:', bot 
!  print*, 'Process:', rank ,' has left:',lft, 'and right:', rgt 
  

! root process spits out
  if(rank .EQ. root) then                                       
     print*, "Program started"
     print*, 'C: ', C
     print*, 'omega: ', omega
     print*, 'trap (ly,lz): ', lambday, lambdaz
     print*, 'init type: ', instate
  endif
 
! count the total number of points in each dimension
  dimx = sum(nptsx(:)) - (nregx - 1) 
  dimy = sum(nptsy(:)) - (nregy - 1) 
  dimz = sum(nptsz(:)) - (nregz - 1) 

! root process spits out
  if(rank .EQ. root) then
     print*, "Points X: ", dimx 
     print*, "Points Y: ", dimy 
     print*, "Points Z: ", dimz 
  endif

! start the CPU counter and open file for saving time, <angular momentum> and <linear momentum>
  if(rank .EQ. root) then
     call cpu_time(tbeg)
     open(unit=23,file= outdir // 'TLP_' // outfile // '.txt',status='replace')
  endif


  !==============================================
  !Define border mesh for parallel communications
  !==============================================
  if(mod(nregx,dims(2)) .NE. 0 .OR. mod(nregy,dims(1)).NE.0) then
     print*, 'Number of dimensions and processes do not match. Quitting.'
     stop
  endif

! set up the indexing vectors
  !==== ROWS:points ===
  x_begind(1) = 1
  do dumind = 1,dims(2)
     x_endind(dumind) = x_begind(dumind) + ceiling(real(dimx) / real(dims(2))) - 1
     if(dumind .LT. dims(2)) then
        x_begind(dumind + 1) = x_endind(dumind) 
     endif
  enddo
  x_endind(dims(2)) = dimx
  !==== ROWS:blocks ===
  x_begblk(1) = 1
  do dumind = 1,dims(2)
     x_endblk(dumind) = x_begblk(dumind) + nregx / dims(2) - 1
     if(dumind .LT. dims(2)) then
        x_begblk(dumind + 1) = x_endblk(dumind) + 1
     endif
  enddo
  x_endblk(dims(2)) = nregx

  !==== COLUMNS:points ===
  y_begind(1) = 1
  do dumind = 1,dims(1)
     y_endind(dumind) = y_begind(dumind) + ceiling(real(dimy) / real(dims(1))) - 1
     if(dumind .LT. dims(1)) then
        y_begind(dumind + 1) = y_endind(dumind) 
     endif
  enddo
  y_endind(dims(1)) = dimy

  !==== COLUMNGS:blocks ===
  y_begblk(1) = 1
  do dumind = 1,dims(1)
     y_endblk(dumind) = y_begblk(dumind) + nregy / dims(1) - 1
     if(dumind .LT. dims(1)) then
        y_begblk(dumind + 1) = y_endblk(dumind) + 1
     endif
  enddo
  y_endblk(dims(1)) = nregy
  !=========================================================================

! set up the "self dimensions" of each process 
  my_dimx = x_endind(my_xind) - x_begind(my_xind) + 1
  my_dimy = y_endind(my_yind) - y_begind(my_yind) + 1
 






  !===============================================================
  !Create basis
  !===============================================================
  if(rank .EQ. root) then
     print*, "Constructing FEM-DVR basis..."
  endif
  call kanta 
  call ptswts
  !===============================================================
  !Initialize wave function
  !===============================================================
  if(rank .EQ. root) then
     print*, "Initializing wavefunction..."
  endif
  call init_wf 
  !===============================================================
  !construct kinetic energy block matrices
  !===============================================================
  if(rank .EQ. root) then
     print*, "Constructing kinetic energy block matrices..."
  endif
  call ke_lohkot
  ! ==============================================================
  ! construct kinetic energy block propagators
  ! ==============================================================
  if(rank .EQ. root) then
     print*, "Constructing kinetic energy block propagators..."
  endif
  call ke_props
  ! ==============================================================
  ! construct momentum block matrices
  ! ==============================================================
  if(rank .EQ. root) then
     print*, "Constructing momentum operator block matrices..."
  endif
  call mom_lohkot
  ! ==============================================================
  ! construct kinetic energy block propagators
  ! ==============================================================
  if(rank .EQ. root) then
     print*, "Constructing momentum operator block propagators..."
  endif
  call mom_props
  ! ==============================================================
  ! start propagating
  ! ============================================================== 
  call virhenormi
  if(rank .EQ. root) then
     print*, "Propagating..."
  endif
  call propagate
  
  if(rank .EQ. root) close(unit=23)
  if(rank .EQ. root) call cpu_time(tend)
  if(rank .EQ. root) print*, 'Propagation finished in ', tend - tbeg ,'CPU seconds.' 
  ! ============ finish MPI =======================================
  call MPI_Finalize(ierr)              

end program GPfemdvr

