
!==============================================================================
! P propagators
!==============================================================================

subroutine p_koddpx
  ! ====================================
  ! odd kinetic energy block propagation
  ! ====================================
  use globaali, only : nregx,nregy,nregz, nptsx,nptsy,nptsz, wf, lohkox,lohkoy,lohkoz
  use globaali, only : dimy, dimz
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, i,j,k
  
  ! first dimension 
  do j = 1, my_dimy
     do k = 1, dimz
        beg = 1 
        fin = beg + nptsx( x_begblk(my_xind) ) - 1  
        do ind = x_begblk(my_xind), x_endblk(my_xind) - 1, 2
           wf(beg:fin,j,k) = matmul(lohkox(ind)%kprop, wf(beg:fin,j,k) ) 
           if(ind .LT. nregx - 1) then
              beg = fin + nptsx(ind + 1) - 1
              fin = beg + nptsx(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine p_koddpx



subroutine p_koddpy
  ! ====================================
  ! odd kinetic energy block propagation
  ! ====================================
  use globaali, only : nregx,nregy,nregz, nptsx,nptsy,nptsz, wf, lohkox,lohkoy,lohkoz
  use globaali, only : dimx, dimz
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, i,j,k
  
  
  ! second dimension 
  do i = 1, my_dimx
     do k = 1, dimz
        beg = 1
        fin = beg        + nptsy( y_begblk(my_yind) ) - 1  
        do ind = y_begblk(my_yind), y_endblk(my_yind) - 1, 2
           wf(i,beg:fin,k) = matmul(lohkoy(ind)%kprop, wf(i,beg:fin,k) ) 
           if(ind .LT. nregy - 1) then
              beg = fin + nptsy(ind + 1) - 1
              fin = beg + nptsy(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine p_koddpy




subroutine p_koddpz
  ! ====================================
  ! odd kinetic energy block propagation
  ! ====================================
  use globaali, only : nregx,nregy,nregz, nptsx,nptsy,nptsz, wf, lohkox,lohkoy,lohkoz
  use globaali, only : dimx, dimy
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, i,j,k
  
  ! third dimension 
  do i = 1, my_dimx
     do j = 1, my_dimy
        beg = 1
        fin = nptsz(1)   
        do ind = 1, nregz - 1, 2
           wf(i,j,beg:fin) = matmul(lohkoz(ind)%kprop, wf(i,j,beg:fin) ) 
           if(ind .LT. nregz - 1) then
              beg = fin + nptsz(ind + 1) - 1
              fin = beg + nptsz(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
  
end subroutine p_koddpz





subroutine p_kevenpx
  
  ! ====================================
  ! even kinetic energy block propagation
  ! ====================================
  use globaali, only : nregx,nregy,nregz, nptsx,nptsy,nptsz, wf, lohkox,lohkoy,lohkoz
  use globaali, only : dimy, dimz
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, i,j,k
  
  ! first dimension 
  do j = 1, my_dimy
     do k = 1, dimz
        beg = nptsx(x_begblk(my_xind)) 
        fin = beg               + nptsx(x_begblk(my_xind) + 1) - 1   
        do ind = x_begblk(my_xind) + 1, x_endblk(my_xind), 2
           wf(beg:fin,j,k) = matmul(lohkox(ind)%kprop, wf(beg:fin,j,k) ) 
           if(ind .LT. nregx - 1) then
              beg = fin + nptsx(ind + 1) - 1
              fin = beg + nptsx(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
  
end subroutine p_kevenpx

subroutine p_kevenpy
  
  ! ====================================
  ! even kinetic energy block propagation
  ! ====================================
  use globaali, only : nregx,nregy,nregz, nptsx,nptsy,nptsz, wf, lohkox,lohkoy,lohkoz
  use globaali, only : dimx, dimz
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, i,j,k
  
  ! second dimension 
  do i = 1, my_dimx
     do k = 1, dimz
        beg = nptsy(y_begblk(my_yind)) 
        fin = beg               + nptsy(y_begblk(my_yind) + 1 ) - 1   
        do ind = y_begblk(my_yind) + 1, y_endblk(my_yind), 2
           wf(i,beg:fin,k) = matmul(lohkoy(ind)%kprop, wf(i,beg:fin,k) ) 
           if(ind .LT. nregy - 1) then
              beg = fin + nptsy(ind + 1) - 1
              fin = beg + nptsy(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine p_kevenpy
  

subroutine p_kevenpz
  
  ! ====================================
  ! even kinetic energy block propagation
  ! ====================================
  use globaali, only : nregx,nregy,nregz, nptsx,nptsy,nptsz, wf, lohkox,lohkoy,lohkoz
  use globaali, only : dimx, dimy
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, i,j,k
  
  ! third dimension 
  do i = 1, my_dimx
     do j = 1, my_dimy
        beg = nptsz(1)
        fin = beg + nptsz(2) - 1   
        do ind = 2, nregz, 2
           wf(i,j,beg:fin) = matmul(lohkoz(ind)%kprop, wf(i,j,beg:fin) ) 
           if(ind .LT. nregz - 1) then
              beg = fin + nptsz(ind + 1) - 1
              fin = beg + nptsz(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
  
end subroutine p_kevenpz



!==============================================================================
! Q propagators
!==============================================================================


subroutine p_koddqx
  ! ====================================
  ! odd kinetic energy block propagation
  ! ====================================
  use globaali, only : nregx,nregy,nregz, nptsx,nptsy,nptsz, wf, lohkox,lohkoy,lohkoz
  use globaali, only : dimy, dimz
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, i,j,k
  
  ! first dimension 
  do j = 1, my_dimy
     do k = 1, dimz
        beg = 1
        fin = beg + nptsx( x_begblk(my_xind) ) - 1  
        do ind = x_begblk(my_xind), x_endblk(my_xind) - 1, 2
           wf(beg:fin,j,k) = matmul(lohkox(ind)%kpropq, wf(beg:fin,j,k) ) 
           if(ind .LT. nregx - 1) then
              beg = fin + nptsx(ind + 1) - 1
              fin = beg + nptsx(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine p_koddqx



subroutine p_koddqy
  ! ====================================
  ! odd kinetic energy block propagation
  ! ====================================
  use globaali, only: nregx,nregy,nregz, nptsx,nptsy,nptsz, wf, lohkox,lohkoy,lohkoz
  use globaali, only: dimx, dimz
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###

  implicit none
  integer :: ind, beg, fin, i,j,k
  
  
  ! second dimension 
  do i = 1, my_dimx
     do k = 1, dimz
        beg = 1 
        fin = beg        + nptsy( y_begblk(my_yind) ) - 1  
        do ind = y_begblk(my_yind), y_endblk(my_yind) - 1, 2
           wf(i,beg:fin,k) = matmul(lohkoy(ind)%kpropq, wf(i,beg:fin,k) ) 
           if(ind .LT. nregy - 1) then
              beg = fin + nptsy(ind + 1) - 1
              fin = beg + nptsy(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine p_koddqy




subroutine p_koddqz
  ! ====================================
  ! odd kinetic energy block propagation
  ! ====================================
  use globaali, only: nregx,nregy,nregz, nptsx,nptsy,nptsz, wf, lohkox,lohkoy,lohkoz
  use globaali, only: dimx, dimy
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, i,j,k
  
  ! third dimension 
  do i = 1, my_dimx
     do j = 1, my_dimy
        beg = 1
        fin = nptsz(1)   
        do ind = 1, nregz - 1, 2
           wf(i,j,beg:fin) = matmul(lohkoz(ind)%kpropq, wf(i,j,beg:fin) ) 
           if(ind .LT. nregz - 1) then
              beg = fin + nptsz(ind + 1) - 1
              fin = beg + nptsz(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
  
end subroutine p_koddqz





subroutine p_kevenqx
  
  ! ====================================
  ! even kinetic energy block propagation
  ! ====================================
  use globaali, only : nregx,nregy,nregz, nptsx,nptsy,nptsz, wf, lohkox,lohkoy,lohkoz
  use globaali, only : dimy, dimz
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, i,j,k
  
  ! first dimension 
  do j = 1, my_dimy
     do k = 1, dimz
        beg = nptsx(x_begblk(my_xind)) 
        fin = beg               + nptsx(x_begblk(my_xind) + 1) - 1   
        do ind = x_begblk(my_xind) + 1, x_endblk(my_xind), 2 
           wf(beg:fin,j,k) = matmul(lohkox(ind)%kpropq, wf(beg:fin,j,k) ) 
           if(ind .LT. nregx - 1) then
              beg = fin + nptsx(ind + 1) - 1
              fin = beg + nptsx(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
  
end subroutine p_kevenqx

subroutine p_kevenqy
  
  ! ====================================
  ! even kinetic energy block propagation
  ! ====================================
  use globaali, only : nregx,nregy,nregz, nptsx,nptsy,nptsz, wf, lohkox,lohkoy,lohkoz
  use globaali, only : dimx, dimz
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, i,j,k
  
  ! second dimension 
  do i = 1, my_dimx
     do k = 1, dimz
        beg = nptsy(y_begblk(my_yind)) 
        fin = beg               + nptsy(y_begblk(my_yind) + 1 ) - 1   
        do ind = y_begblk(my_yind) + 1, y_endblk(my_yind), 2 
           wf(i,beg:fin,k) = matmul(lohkoy(ind)%kpropq, wf(i,beg:fin,k) ) 
           if(ind .LT. nregy - 1) then
              beg = fin + nptsy(ind + 1) - 1
              fin = beg + nptsy(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine p_kevenqy
  

subroutine p_kevenqz
  
  ! ====================================
  ! even kinetic energy block propagation
  ! ====================================
  use globaali, only : nregx,nregy,nregz, nptsx,nptsy,nptsz, wf, lohkox,lohkoy,lohkoz
  use globaali, only : dimx, dimy
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  integer :: ind, beg, fin, i,j,k
  
  ! third dimension 
  do i = 1, my_dimx
     do j = 1, my_dimy
        beg = nptsz(1)
        fin = beg + nptsz(2) - 1   
        do ind = 2, nregz, 2
           wf(i,j,beg:fin) = matmul(lohkoz(ind)%kpropq, wf(i,j,beg:fin) ) 
           if(ind .LT. nregz - 1) then
              beg = fin + nptsz(ind + 1) - 1
              fin = beg + nptsz(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
  
end subroutine p_kevenqz
