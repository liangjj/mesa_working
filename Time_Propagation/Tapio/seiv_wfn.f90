subroutine seiv_wfn(nro)
  !===================================================
  ! save wavefunction
  ! ==================================================
  use globaali, only: io, outdir, outfile, pwx, pwy, pwz, wf
  use globaali, only: dimx, dimy, dimz
  use globaali, only: rank, root, ierr, size                            ! ### MPI ###
  use globaali, only: cart_comm, dims                                   ! ### MPI ###
  use globaali, only: top, bot, lft, rgt                                ! ### MPI ###
  use globaali, only: x_begind, x_endind, y_begind, y_endind            ! ### MPI ###
  use globaali, only: my_xind, my_yind                                  ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                                  ! ### MPI ###
  implicit none
  include 'mpif.h'

  character(LEN=3) :: cnro
  integer          :: ones, tens, dogs
  integer          :: status(MPI_Status_Size), i_xind, i_yind
  integer          :: nro, i, j, k, count, countx, county

  complex*16,allocatable, dimension(:,:,:) :: savewf
  if(rank.EQ.root) then
     allocate(savewf(dimx,dimy,dimz))
     savewf = 0d0
     savewf( x_begind(my_xind):x_endind(my_xind), y_begind(my_yind):y_endind(my_yind),:) = wf
  endif

  countx = my_dimx 
  county = my_dimy
  count  = countx * county * dimz
  
  call MPI_Barrier(cart_comm,ierr)
  !========== gather parts of wavefunction from processes to root for saving =====================
  if(rank .NE. root) then
     call MPI_Send(wf, count, MPI_Double_Complex, root, rank, cart_comm, ierr)
  else if(rank.EQ.root) then
     do i = 1,size-1
        i_xind = 1 + mod(i,dims(2))
        i_yind = ceiling((i+1) /real(dims(2)))
        countx = (x_endind(i_xind) - x_begind(i_xind)) + 1
        county = (y_endind(i_yind) - y_begind(i_yind)) + 1
        count  = countx * county * dimz
        call MPI_Recv(savewf( x_begind(i_xind):x_endind(i_xind), y_begind(i_yind):y_endind(i_yind),:), count, MPI_Double_Complex, i, i,  cart_comm, status, ierr)
     enddo
  endif
  

  !================= SAVE =======================
  if (io .AND. rank .EQ. root) then
     ones = mod(nro,10)
     tens = mod(nro - ones,100)/10
     dogs = mod(nro - ones - tens, 1000)/100
     ones = ones + 48
     tens = tens + 48
     dogs = dogs + 48
     cnro = char(dogs) // char(tens) // char(ones)


     open(unit=11,file= outdir // outfile // '_' // cnro //'.wfn')
     do i = 1, dimx
        do j = 1, dimy
           do k = 1, dimz
              write(11,100) real(savewf(i,j,k)), aimag(savewf(i,j,k))
           enddo
        enddo
     enddo
     close(unit=11)
  end if
100 format (E13.5E3,1X,E13.5E3)


if(rank.EQ.root) deallocate(savewf)
end subroutine seiv_wfn
