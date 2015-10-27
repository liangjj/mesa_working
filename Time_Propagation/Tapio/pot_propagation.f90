
subroutine p_vdiag
  ! =========================================
  ! apply diagonal (potential) propagator
  ! =========================================
  use globaali, only: wf, vprop
  use globaali, only: dimx,dimy,dimz
  use globaali, only: x_begind, x_endind, y_begind, y_endind            ! ### MPI ###
  use globaali, only: my_xind, my_yind                                  ! ### MPI ###
  use globaali, only : my_dimx, my_dimy                                 ! ### MPI ###
  implicit none
  integer :: i,j,k 
  call ve_prop
  do i = 1,my_dimx
     do j = 1, my_dimy 
        do k = 1, dimz
           wf(i,j,k) = vprop(i,j,k)*wf(i,j,k)
        enddo
     enddo
  enddo
end subroutine p_vdiag




subroutine ve_prop
! =========================================
! update potential and potential propagator
! =========================================
  use globaali, only: C, pq, aika, dt, iu, io, pwx,pwy,pwz, wf, vext, vprop, lambdax, lambday, lambdaz, gamma,myy
  use globaali, only: dimx, dimy, dimz
  use globaali, only: x_begind, x_endind, y_begind, y_endind                  ! ### MPI ###
  use globaali, only: my_xind, my_yind, rank, cart_comm,ierr                  ! ### MPI ###
  implicit none	    
  
  integer  :: i,j,k,myi, myj
 
  vext  = 0d0
  vprop = 0d0
  myi   = 0
  myj   = 0
  do i = x_begind(my_xind), x_endind(my_xind) 
     myi = myi + 1
     myj = 0 
     do j = y_begind(my_yind), y_endind(my_yind)
        myj = myj + 1
        do k = 1, dimz
           ! =====================================   harmonic trap  ================================================
           vext(myi,myj,k)   = 0.5d0 * ((lambdax*pwx(i,1))**2 + (lambday * pwy(j,1))**2 + (lambdaz * pwz(k,1))**2 )         
           ! =====================================   |psi|^2  ======================================================
           vext(myi,myj,k)   = vext(myi,myj,k) + C * abs( wf(myi,myj,k) / sqrt( pwx(i,2)*pwy(j,2)*pwz(k,2)) )**2d0    
           ! =====================================   add here your own time dependent potentials ===================
           !
           !           
           ! =======================================================================================================
          

           ! ==========================   the propagator itself ====================================================
           vprop(myi,myj,k)  = exp( 1/(iu-gamma) * pq(3) * dt * (vext(myi,myj,k) + iu*gamma*myy ) / 2d0 )                 

        enddo
     enddo
  enddo

end subroutine ve_prop

