subroutine normalize
  ! =========================================
  ! normalizes wavefunction
  ! =========================================
  use globaali, only: wf
  use globaali, only: dimx, dimy, dimz
  use globaali, only: root, rank, ierr                       ! ### MPI ###
  use globaali, only: cart_comm                              ! ### MPI ###
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: bot, rgt                               ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none	    
  include 'mpif.h'
  
  integer  :: i, j, k, xend, yend
  real*8   :: normi, totalnorm


  if(bot .NE. MPI_Proc_Null) yend = my_dimy - 1
  if(bot .EQ. MPI_Proc_Null) yend = my_dimy
  if(rgt .NE. MPI_Proc_Null) xend = my_dimx - 1
  if(rgt .EQ. MPI_Proc_Null) xend = my_dimx

  !==== compute local norm =========
  normi = 0d0
  do i = 1, xend
     do j = 1, yend
        do k = 1, dimz
           normi = normi + abs(wf(i,j,k))**2
        enddo
     enddo
  enddo


  !==== root computes total norm ==========
  call MPI_Reduce(normi, totalnorm, 1, MPI_Double_Precision, MPI_Sum, root, cart_comm, ierr)
  !==== root distributes the result ======= 
  call MPI_Bcast(totalnorm, 1, MPI_Double_Precision, root, cart_comm, ierr)
  !==== normalize =========================

  wf = wf / sqrt(totalnorm)

end subroutine normalize

