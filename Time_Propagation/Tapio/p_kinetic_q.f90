subroutine p_kinetic_q
  use globaali, only: wf 
  use globaali, only: dimz
  use globaali, only: rank, ierr                                        ! ### MPI ###
  use globaali, only: cart_comm                                         ! ### MPI ###
  use globaali, only: top, bot, lft, rgt                                ! ### MPI ###
  use globaali, only: x_begind, x_endind, y_begind, y_endind            ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk            ! ### MPI ###
  use globaali, only: my_xind, my_yind                                  ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                                  ! ### MPI ###
  implicit none
  include 'mpif.h' 

  integer :: status(MPI_Status_Size)
  integer :: Ax, Bx, Ay, By
  integer :: count, countx, county
  Ax     = 1
  Bx     = my_dimx
  Ay     = 1
  By     = my_dimy
  countx = Bx 
  county = By 

  call p_koddqx
  !=== send to left ===
  call MPI_Send(wf(Ax,Ay:By,:), county*dimz, MPI_Double_Complex, lft, rank, cart_comm, ierr)
  call MPI_Recv(wf(Bx,Ay:By,:), county*dimz, MPI_Double_Complex, rgt, rgt,  cart_comm, status, ierr)
  !====================
  call p_kevenqx
  !=== send to right ==
  call MPI_Send(wf(Bx,Ay:By,:), county*dimz, MPI_Double_Complex, rgt, rank, cart_comm, ierr)
  call MPI_Recv(wf(Ax,Ay:By,:), county*dimz, MPI_Double_Complex, lft, lft,  cart_comm, status, ierr)
  !====================
  call p_koddqx
  !=== send to left ===
  call MPI_Send(wf(Ax,Ay:By,:), county*dimz, MPI_Double_Complex, lft, rank, cart_comm, ierr)
  call MPI_Recv(wf(Bx,Ay:By,:), county*dimz, MPI_Double_Complex, rgt, rgt,  cart_comm, status, ierr)
  !====================
  call p_koddqy
  !=== send to up =====
  call MPI_Send(wf(Ax:Bx,Ay,:), countx*dimz, MPI_Double_Complex, top, rank, cart_comm, ierr)
  call MPI_Recv(wf(Ax:Bx,By,:), countx*dimz, MPI_Double_Complex, bot, bot,  cart_comm, status, ierr)
  !====================
  call p_kevenqy
  !=== send to down ===
  call MPI_Send(wf(Ax:Bx,By,:), countx*dimz, MPI_Double_Complex, bot, rank, cart_comm, ierr)
  call MPI_Recv(wf(Ax:Bx,Ay,:), countx*dimz, MPI_Double_Complex, top, top,  cart_comm, status, ierr)
  !====================
  call p_koddqy
  !=== send to up =====
  call MPI_Send(wf(Ax:Bx,Ay,:), countx*dimz, MPI_Double_Complex, top, rank, cart_comm, ierr)
  call MPI_Recv(wf(Ax:Bx,By,:), countx*dimz, MPI_Double_Complex, bot, bot,  cart_comm, status, ierr)
  !====================
  
  call p_koddqz
  call p_kevenqz
  call p_koddqz

end subroutine p_kinetic_q
