! ====================================
! P propagation
! ====================================

subroutine p_xdy_odd
  ! ====================================
  ! xdy propagation
  ! ====================================
  use globaali, only : nregy, nptsy, wf, lyprop
  use globaali, only: dimx, dimz
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, i, k, myi
  ! along x dimension 
  do i = 1, my_dimx
     myi = x_begind(my_xind)- 1 + i
     ! odd x*dy block
     do k = 1, dimz
        beg = 1 
        fin = beg        + nptsy( y_begblk(my_yind) ) - 1  
        do ind = y_begblk(my_yind), y_endblk(my_yind) - 1, 2
           wf(i,beg:fin,k) = matmul(lyprop(myi)%blk(ind)%mprop, wf(i,beg:fin,k) ) 
           if(ind .LT. nregy - 1) then
              beg = fin + nptsy(ind + 1) - 1
              fin = beg + nptsy(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine p_xdy_odd
   

subroutine p_xdy_even
  ! ====================================
  ! xdy propagation
  ! ====================================
  use globaali, only : nregy, nptsy, wf, lyprop
  use globaali, only: dimx, dimz
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, i, k, myi
  ! along x dimension 
  do i = 1, my_dimx
     myi = x_begind(my_xind)- 1 + i
     ! even x*dy block
     do k = 1, dimz
        beg = nptsy(y_begblk(my_yind)) 
        fin = beg               + nptsy(y_begblk(my_yind) + 1 ) - 1   
        do ind = y_begblk(my_yind) + 1, y_endblk(my_yind), 2
           wf(i,beg:fin,k) = matmul(lyprop(myi)%blk(ind)%mprop, wf(i,beg:fin,k) ) 
           if(ind .LT. nregy - 1) then
              beg = fin + nptsy(ind + 1) - 1
              fin = beg + nptsy(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine p_xdy_even


subroutine p_ydx_odd
  ! ====================================
  ! ydx propagation
  ! ====================================
  use globaali, only : nregx, nptsx, wf, lxprop
  use globaali, only: dimy, dimz
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, j,k, myj
  ! along y dimension 
  do j = 1, my_dimy
     myj = y_begind(my_yind)- 1 + j
     ! odd y*dx block 
     do k = 1, dimz
        beg = 1 
        fin = beg + nptsx( x_begblk(my_xind) ) - 1  
        do ind = x_begblk(my_xind), x_endblk(my_xind) - 1, 2
           wf(beg:fin,j,k) = matmul(lxprop(myj)%blk(ind)%mprop, wf(beg:fin,j,k) ) 
           if(ind .LT. nregx - 1) then
              beg = fin + nptsx(ind + 1) - 1
              fin = beg + nptsx(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine p_ydx_odd



subroutine p_ydx_even
  ! ====================================
  ! ydx propagation
  ! ====================================
  use globaali, only: nregx, nptsx, wf, lxprop
  use globaali, only: dimy, dimz
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, j,k, myj
  ! along y dimension 
  do j = 1, my_dimy
     myj = y_begind(my_yind)- 1 + j
     ! even y*dx block 
     do k = 1, dimz
        beg = nptsx(x_begblk(my_xind)) 
        fin = beg               + nptsx(x_begblk(my_xind) + 1) - 1   
        do ind = x_begblk(my_xind) + 1, x_endblk(my_xind), 2
           wf(beg:fin,j,k) = matmul(lxprop(myj)%blk(ind)%mprop, wf(beg:fin,j,k) ) 
           if(ind .LT. nregx - 1) then
              beg = fin + nptsx(ind + 1) - 1
              fin = beg + nptsx(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine p_ydx_even


! ==============================================================================
! =========== Q propagation ====================================================
! ==============================================================================


subroutine p_xdy_oddq
  ! ====================================
  ! xdy propagation
  ! ====================================
  use globaali, only : nregy, nptsy, wf, lyprop
  use globaali, only: dimx, dimz
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, i, k,myi
  ! along x dimension 
  do i = 1, my_dimx
     myi = x_begind(my_xind)- 1 + i
     ! odd x*dy block
     do k = 1, dimz
        beg = 1 
        fin = beg        + nptsy( y_begblk(my_yind) ) - 1  
        do ind = y_begblk(my_yind), y_endblk(my_yind) - 1, 2
           wf(i,beg:fin,k) = matmul(lyprop(myi)%blk(ind)%mpropq, wf(i,beg:fin,k) ) 
           if(ind .LT. nregy - 1) then
              beg = fin + nptsy(ind + 1) - 1
              fin = beg + nptsy(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine p_xdy_oddq
   

subroutine p_xdy_evenq
  ! ====================================
  ! xdy propagation
  ! ====================================
  use globaali, only : nregy, nptsy, wf, lyprop
  use globaali, only: dimx, dimz
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, i, k, myi
  ! along x dimension 
  do i = 1, my_dimx
     myi = x_begind(my_xind)- 1 + i
     ! even x*dy block
     do k = 1, dimz
        beg = nptsy(y_begblk(my_yind)) 
        fin = beg               + nptsy(y_begblk(my_yind) + 1 ) - 1   
        do ind = y_begblk(my_yind) + 1, y_endblk(my_yind), 2
           wf(i,beg:fin,k) = matmul(lyprop(myi)%blk(ind)%mpropq, wf(i,beg:fin,k) ) 
           if(ind .LT. nregy - 1) then
              beg = fin + nptsy(ind + 1) - 1
              fin = beg + nptsy(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine p_xdy_evenq


subroutine p_ydx_oddq
  ! ====================================
  ! ydx propagation
  ! ====================================
  use globaali, only : nregx, nptsx, wf, lxprop
  use globaali, only: dimy, dimz
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, j,k, myj
  ! along y dimension 
  do j = 1, my_dimy
     myj = y_begind(my_yind)- 1 + j
     ! odd y*dx block 
     do k = 1, dimz
        beg = 1 
        fin = beg + nptsx( x_begblk(my_xind) ) - 1  
        do ind = x_begblk(my_xind), x_endblk(my_xind) - 1, 2
           wf(beg:fin,j,k) = matmul(lxprop(myj)%blk(ind)%mpropq, wf(beg:fin,j,k) ) 
           if(ind .LT. nregx - 1) then
              beg = fin + nptsx(ind + 1) - 1
              fin = beg + nptsx(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine p_ydx_oddq



subroutine p_ydx_evenq
  ! ====================================
  ! ydx propagation
  ! ====================================
  use globaali, only: nregx, nptsx, wf, lxprop
  use globaali, only: dimy, dimz
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, j,k, myj
  ! along y dimension 
  do j = 1, my_dimy
     myj = y_begind(my_yind)- 1 + j
     ! even y*dx block 
     do k = 1, dimz
        beg = nptsx(x_begblk(my_xind)    ) 
        fin = beg               + nptsx(x_begblk(my_xind) + 1) - 1   
        do ind = x_begblk(my_xind) + 1, x_endblk(my_xind), 2
           wf(beg:fin,j,k) = matmul(lxprop(myj)%blk(ind)%mpropq, wf(beg:fin,j,k) ) 
           if(ind .LT. nregx - 1) then
              beg = fin + nptsx(ind + 1) - 1
              fin = beg + nptsx(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine p_ydx_evenq

