! ==================================================================
! This file contains subroutines for computing various expectation
! values from the wavefunction. Also the error norm is computed here
! ==================================================================
subroutine virhenormi
  ! =========================================
  ! update potential and potential propagator
  ! =========================================
  use globaali, only: outdir, outfile
  use globaali, only: C, omega, myy, ang, err, aika, lambdax, lambday, lambdaz, iu
  use globaali, only: wf, hwf, hpsi, vext, locerr, ehwf
  use globaali, only: pwx, pwy, pwz, dimx, dimy, dimz
  use globaali, only: root, rank, ierr, size                 ! ### MPI ###
  use globaali, only: cart_comm, dims                        ! ### MPI ###
  use globaali, only: top, bot, lft, rgt                     ! ### MPI ###
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###

  implicit none	    
  include 'mpif.h'  
  integer         :: status(MPI_Status_Size)

  logical         :: fileflag
  integer         :: i, j, k, mind, mmax
  real*8          :: pi, kin, kod, kev, pot, int, ploc, intloc, virhe, fac
  real*8          :: totalkod, totalkev, totalpot, totalint, err_rank
  complex*16      :: pod, pev, totalpod, totalpev
  complex*16      :: ydx, xdy, totalydx, totalxdy
  integer         :: Ax, Bx, Ay, By, xend, yend, countx, county, count, i_xind,i_yind
  integer         :: Aix, Bix, Aiy, Biy, myi, myj
  complex*16      :: linmom, cmloca, cmlocb,cmtot
  complex*16 ,allocatable,dimension(:,:) :: cma,cmb
  



  ! ALL BLOCKS NEED TO BE SAME SIZE!!!

  call MPI_Barrier(cart_comm,ierr)

  Ax = 1
  Bx = my_dimx
  Ay = 1
  By = my_dimy
  if(bot .NE. MPI_Proc_Null) yend = By - 1
  if(bot .EQ. MPI_Proc_Null) yend = By
  if(rgt .NE. MPI_Proc_Null) xend = Bx - 1
  if(rgt .EQ. MPI_Proc_Null) xend = Bx
  
  kod = 0d0
  kev = 0d0
  kin = 0d0
  pot = 0d0
  int = 0d0
  xdy = 0d0
  ydx = 0d0
  ang = 0d0
  pod = 0d0
  pev = 0d0
  linmom = 0d0
  totalpod = 0d0 
  totalpev = 0d0
  totalpot = 0d0
  totalkev = 0d0
  totalkod = 0d0
  totalint = 0d0
  err_rank = 0d0
  err      = 0d0
  
  ! ================================================================================
  ! ================================ kinetic =======================================
  ! ================================================================================
  hpsi  = 0d0
  hwf   = 0d0
  ehwf  = 0d0
  call koddx_wf
  call koddy_wf
  call koddz_wf
  hpsi = hpsi + ehwf
  do i = Ax, Bx 
    do j = Ay, By
       do k = 1, dimz        
          kod = kod + hwf(i,j,k) * conjg(wf(i,j,k)) 
      enddo
     enddo
  enddo
  !==== root computes total sum ==========
  call MPI_Reduce(kod, totalkod, 1, MPI_Double_Precision, MPI_Sum, root, cart_comm, ierr)
  
  hwf   = 0d0
  ehwf  = 0d0
  call kevenx_wf
  call keveny_wf
  call kevenz_wf
  hpsi = hpsi + ehwf
  do i = Ax, Bx 
     do j = Ay, By 
        do k = 1, dimz
           kev = kev + hwf(i,j,k) * conjg(wf(i,j,k))            
       enddo
     enddo
  enddo
  !==== root computes total sum ==========
  call MPI_Reduce(kev, totalkev, 1, MPI_Double_Precision, MPI_Sum, root, cart_comm, ierr)


if(1.EQ.1) then
   hwf = 0d0
  ! ================================================================================
  ! ================================ P_z  =======================================
  ! ================================================================================
  call oddpz_wf
   do i = Ax, Bx 
    do j = Ay, By
       do k = 1, dimz        
          pod = pod + hwf(i,j,k) * conjg(wf(i,j,k))  
      enddo
     enddo
  enddo
  !==== root computes total sum ==========
  call MPI_Reduce(pod, totalpod, 1, MPI_Double_Complex, MPI_Sum, root, cart_comm, ierr)
  hwf =0d0
  call evenpz_wf
  do i = Ax, Bx 
     do j = Ay, By 
        do k = 1, dimz
           pev = pev + hwf(i,j,k) * conjg(wf(i,j,k))          
       enddo
     enddo
  enddo
  !==== root computes total sum ==========
  call MPI_Reduce(pev, totalpev, 1, MPI_Double_Complex, MPI_Sum, root, cart_comm, ierr)
endif
linmom = -iu*(totalpod + totalpev)


  ! ===================================================================================
  ! =================================== potential ===================================== 
  ! ===================================================================================  
  do i = Ax, Bx
     do j = Ay, By
        do k = 1, dimz
           fac = 1d0
           if(i.EQ.1 .AND.lft.NE.MPI_Proc_Null) fac = 1/2d0
           if(i.EQ.Bx.AND.rgt.NE.MPI_Proc_Null) fac = 1/2d0
           if(j.EQ.1 .AND.top.NE.MPI_Proc_Null) fac = 1/2d0
           if(j.EQ.By.AND.bot.NE.MPI_Proc_Null) fac = 1/2d0

           if(i.EQ.1 .AND. j.EQ.1   .AND. lft.NE.MPI_Proc_Null .AND. top .NE.MPI_Proc_Null) fac = 1/4d0
           if(i.EQ.1 .AND. j.EQ.By  .AND. lft.NE.MPI_Proc_Null .AND. bot .NE.MPI_Proc_Null) fac = 1/4d0
           if(i.EQ.Bx .AND. j.EQ.1  .AND. rgt.NE.MPI_Proc_Null .AND. top .NE.MPI_Proc_Null) fac = 1/4d0
           if(i.EQ.Bx .AND. j.EQ.By .AND. rgt.NE.MPI_Proc_Null .AND. bot .NE.MPI_Proc_Null) fac = 1/4d0

           ploc = 0.5d0* ((lambdax*pwx(x_begind(my_xind)-1+i,1))**2 + (lambday*pwy(y_begind(my_yind)-1+j,1))**2 + (lambdaz*pwz(k,1))**2  ) 
           pot = pot + ploc *  abs(wf(i,j,k))**2 * fac
           hpsi(i,j,k) = hpsi(i,j,k) + ploc * wf(i,j,k) 
        enddo
     enddo
  enddo
  !==== root computes total sum ==========
  call MPI_Reduce(pot, totalpot, 1, MPI_Double_Precision, MPI_Sum, root, cart_comm, ierr)
  
  ! ====================================================================================== 
  ! ==============================  interaction ==========================================
  ! ====================================================================================== 
  do i = Ax, Bx
     do j = Ay, By
        do k = 1, dimz
           fac = 1d0 
           if(i.EQ.1 .AND.lft.NE.MPI_Proc_Null) fac = 1/2d0
           if(i.EQ.Bx.AND.rgt.NE.MPI_Proc_Null) fac = 1/2d0
           if(j.EQ.1 .AND.top.NE.MPI_Proc_Null) fac = 1/2d0
           if(j.EQ.By.AND.bot.NE.MPI_Proc_Null) fac = 1/2d0

           if(i.EQ.1 .AND. j.EQ.1   .AND. lft.NE.MPI_Proc_Null .AND. top .NE.MPI_Proc_Null) fac = 1/4d0
           if(i.EQ.1 .AND. j.EQ.By  .AND. lft.NE.MPI_Proc_Null .AND. bot .NE.MPI_Proc_Null) fac = 1/4d0
           if(i.EQ.Bx .AND. j.EQ.1  .AND. rgt.NE.MPI_Proc_Null .AND. top .NE.MPI_Proc_Null) fac = 1/4d0
           if(i.EQ.Bx .AND. j.EQ.By .AND. rgt.NE.MPI_Proc_Null .AND. bot .NE.MPI_Proc_Null) fac = 1/4d0
           
           intloc = 0.5d0 * C * abs(wf(i,j,k) / sqrt(pwx(x_begind(my_xind)-1+i,2)*pwy(y_begind(my_yind)-1+j,2)*pwz(k,2)))**2 
           int = int + intloc * abs(wf(i,j,k))**2 * fac
           hpsi(i,j,k) = hpsi(i,j,k) + 2d0 * intloc * wf(i,j,k) 
        enddo
     enddo
  enddo
  !==== root computes total sum ==========
  call MPI_Reduce(int, totalint, 1, MPI_Double_Precision, MPI_Sum, root, cart_comm, ierr)
  

  ! =====================================================================================
  ! =================================== angular =========================================
  ! =====================================================================================
  ehwf = 0d0
  hwf = 0d0
  xdy = 0d0
  call xdy_wf
  hpsi = hpsi - omega * (-iu)* ehwf
  do i = Ax, Bx 
     do j = Ay, By 
        do k = 1, dimz
           xdy = xdy + hwf(i,j,k) * conjg(wf(i,j,k)) 
        enddo
     enddo
  enddo
  !==== root computes total sum ==========
  call MPI_Reduce(xdy, totalxdy, 1, MPI_Double_Complex, MPI_Sum, root, cart_comm, ierr)

  ehwf = 0d0
  hwf = 0d0
  ydx = 0d0
  call ydx_wf
  hpsi = hpsi - omega* (iu) * ehwf
  do i = Ax, Bx 
     do j = Ay, By 
        do k = 1, dimz
           ydx = ydx + hwf(i,j,k) * conjg(wf(i,j,k)) 
        enddo
     enddo
  enddo
 !==== root computes total sum ==========
 call MPI_Reduce(ydx, totalydx, 1, MPI_Double_Complex, MPI_Sum, root, cart_comm, ierr)



 err   = 0d0
 virhe = 0d0
 !==== ERROR ==========
  do i = Ax, Bx 
     do j = Ay, By 
        do k = 1, dimz
           virhe = virhe + abs(hpsi(i,j,k) - myy * wf(i,j,k))**2
       enddo
     enddo
  enddo
!print*,'Virhe, rank: ', virhe, rank
call MPI_Reduce(virhe, err, 1, MPI_Double_Precision, MPI_Sum, root, cart_comm, ierr)

 

  if(rank .EQ. root) then
     kin = totalkod + totalkev
     ang = - iu * (totalxdy - totalydx)
     if(aimag(ang) .GT. 1e-6) then
        print*, 'L_z has imaginary part: ', aimag(ang)
     end if
     
     myy = kin + totalpot + 2d0 * totalint - omega * ang
     err = sqrt(err / myy**2)

  endif

  !==== root distributes chemical potential and error norm ======= 
  call MPI_Bcast(myy, 1, MPI_Double_Precision, root, cart_comm, ierr)
  call MPI_Bcast(err, 1, MPI_Double_Precision, root, cart_comm, ierr)


! ===== write T, L_z, and P_z  
  if(rank.EQ. root) then
     inquire(23,EXIST = fileflag)
     if(fileflag) then
        open(unit=23,file= outdir // 'TLP_' // outfile // '.txt', position='append')
     endif
100  format (E12.5,1X,E12.5,1X,E12.5, 1X, E12.5)

     write(23,100) abs(aika), abs(ang), abs(linmom)
     close(23)
  endif


end subroutine virhenormi







subroutine koddx_wf
  ! ====================================
  ! odd kinetic energy block propagation
  ! ====================================
  use globaali, only: nregx, nregy, nregz
  use globaali, only: nptsx, nptsy, nptsz 
  use globaali, only: wf, hwf, ehwf 
  use globaali, only: lohkox, lohkoy, lohkoz
  use globaali, only: dimy, dimz
  use globaali, only: bot                               ! ### MPI ###
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind, rank                 ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###

  implicit none
  include 'mpif.h'

  integer :: ind, beg, fin, i,j,k
  integer :: Ax, Ay, By, yend

  Ax = 1 
  Ay = 1
  By = my_dimy

  
  ! first dimension 

  if(bot .NE. MPI_Proc_Null) yend = By - 1
  if(bot .EQ. MPI_Proc_Null) yend = By
  do j = Ay, By
     do k = 1, dimz
        beg = Ax 
        fin = beg + nptsx( x_begblk(my_xind) ) - 1
        do ind = x_begblk(my_xind), x_endblk(my_xind) - 1, 2
           if(j.LE.yend) hwf(beg:fin,j,k) = hwf(beg:fin,j,k) + matmul(lohkox(ind)%ke, wf(beg:fin,j,k) ) 
           ehwf(beg:fin,j,k) = ehwf(beg:fin,j,k) + matmul(lohkox(ind)%ke, wf(beg:fin,j,k) ) 
           if(ind .LT. nregx - 1) then
              beg = fin + nptsx(ind + 1) - 1
              fin = beg + nptsx(ind + 2) - 1
           endif
        enddo
     enddo
  enddo

end subroutine koddx_wf


subroutine koddy_wf
  ! ====================================
  ! odd kinetic energy block propagation
  ! ====================================
  use globaali, only: nregx, nregy, nregz 
  use globaali, only: nptsx, nptsy, nptsz 
  use globaali, only: wf, hwf, ehwf 
  use globaali, only: lohkox, lohkoy, lohkoz
  use globaali, only: dimx, dimz 
  use globaali, only: rgt                                    ! ### MPI ###
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###

  implicit none
  include 'mpif.h'

  integer :: ind, beg, fin, i,j,k
  integer :: Ax, Bx, Ay, xend

  Ax = 1 
  Bx = my_dimx
  Ay = 1


 ! second dimension 
  if(rgt .NE. MPI_Proc_Null) xend = Bx - 1 
  if(rgt .EQ. MPI_Proc_Null) xend = Bx
  do i = Ax, Bx
     do k = 1, dimz
        beg = Ay 
        fin = beg + nptsy( y_begblk(my_yind) ) - 1  
        do ind = y_begblk(my_yind), y_endblk(my_yind) - 1, 2
           if(i.LE.xend) hwf(i,beg:fin,k) = hwf(i,beg:fin,k) + matmul(lohkoy(ind)%ke, wf(i,beg:fin,k) ) 
           ehwf(i,beg:fin,k) = ehwf(i,beg:fin,k) + matmul(lohkoy(ind)%ke, wf(i,beg:fin,k) ) 
           if(ind .LT. nregy - 1) then
              beg = fin + nptsy(ind + 1) - 1
              fin = beg + nptsy(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine koddy_wf


subroutine koddz_wf
  ! ====================================
  ! odd kinetic energy block propagation
  ! ====================================
  use globaali, only: nregx, nregy, nregz
  use globaali, only: nptsx, nptsy, nptsz
  use globaali, only: wf, hwf, ehwf 
  use globaali, only: lohkox, lohkoy, lohkoz
  use globaali, only: dimx, dimy 
  use globaali, only: bot, rgt                               ! ### MPI ###
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  include 'mpif.h'

  integer :: ind, beg, fin, i,j,k
  integer :: Ax, Bx, Ay, By, xend, yend

  Ax = 1 
  Bx = my_dimx
  Ay = 1
  By = my_dimy


! third dimension 
  if(bot .NE. MPI_Proc_Null) yend = By - 1 
  if(bot .EQ. MPI_Proc_Null) yend = By
  if(rgt .NE. MPI_Proc_Null) xend = Bx - 1 
  if(rgt .EQ. MPI_Proc_Null) xend = Bx
  do i = Ax, Bx
     do j = Ay, By
        beg = 1
        fin = nptsz(1)   
        do ind = 1, nregz - 1, 2
           if(i.LE.xend .AND. j.LE.yend) hwf(i,j,beg:fin) = hwf(i,j,beg:fin) + matmul(lohkoz(ind)%ke, wf(i,j,beg:fin) ) 
           ehwf(i,j,beg:fin) = ehwf(i,j,beg:fin) + matmul(lohkoz(ind)%ke, wf(i,j,beg:fin) ) 
           if(ind .LT. nregz - 1) then
              beg = fin + nptsz(ind + 1) - 1
              fin = beg + nptsz(ind + 2) - 1
           endif
        enddo
     enddo
  enddo

end subroutine koddz_wf





subroutine kevenx_wf
  
  ! ====================================
  ! even kinetic energy block propagation
  ! ====================================
  use globaali, only: nregx, nregy, nregz 
  use globaali, only: nptsx, nptsy, nptsz 
  use globaali, only: wf, hwf, ehwf 
  use globaali, only: lohkox, lohkoy, lohkoz
  use globaali, only: dimy, dimz
  use globaali, only: bot                                    ! ### MPI ###
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  include 'mpif.h'
 
  integer :: ind, beg, fin, i,j,k
  integer :: Ax, Ay, By, yend

  Ax = 1 
  Ay = 1
  By = my_dimy
  

  ! first dimension 
  if(bot .NE. MPI_Proc_Null) yend = By - 1 
  if(bot .EQ. MPI_Proc_Null) yend = By
  do j = Ay, By
     do k = 1, dimz
        beg = Ax  + nptsx(x_begblk(my_xind)    ) - 1
        fin = beg + nptsx(x_begblk(my_xind) + 1) - 1   
        do ind = x_begblk(my_xind) + 1, x_endblk(my_xind), 2
           if(j.LE.yend) hwf(beg:fin,j,k) = hwf(beg:fin,j,k) + matmul(lohkox(ind)%ke, wf(beg:fin,j,k) ) 
           ehwf(beg:fin,j,k) = ehwf(beg:fin,j,k) + matmul(lohkox(ind)%ke, wf(beg:fin,j,k) ) 
           if(ind .LT. nregx - 1) then
              beg = fin + nptsx(ind + 1) - 1
              fin = beg + nptsx(ind + 2) - 1
           endif
        enddo
     enddo
  enddo

end subroutine kevenx_wf

subroutine keveny_wf
  ! ====================================
  ! even kinetic energy block propagation
  ! ====================================
  use globaali, only: nregx, nregy, nregz 
  use globaali, only: nptsx, nptsy, nptsz 
  use globaali, only: wf, hwf,ehwf 
  use globaali, only: lohkox, lohkoy, lohkoz
  use globaali, only: dimx, dimz
  use globaali, only: rgt                                    ! ### MPI ###
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  include 'mpif.h'

  integer :: ind, beg, fin, i,j,k
  integer :: Ax, Bx, Ay, xend

  Ax = 1 
  Bx = my_dimx
  Ay = 1

  ! second dimension 
  if(rgt .NE. MPI_Proc_Null) xend = Bx - 1 
  if(rgt .EQ. MPI_Proc_Null) xend = Bx
  do i = Ax, Bx
     do k = 1, dimz
        beg = Ay + nptsy(y_begblk(my_yind)     ) - 1
        fin = beg               + nptsy(y_begblk(my_yind) + 1 ) - 1   
        do ind = y_begblk(my_yind) + 1, y_endblk(my_yind), 2
           if(i.LE.xend) hwf(i,beg:fin,k) = hwf(i,beg:fin,k) + matmul(lohkoy(ind)%ke, wf(i,beg:fin,k) ) 
           ehwf(i,beg:fin,k) = ehwf(i,beg:fin,k) + matmul(lohkoy(ind)%ke, wf(i,beg:fin,k) ) 
           if(ind .LT. nregy - 1) then
              beg = fin + nptsy(ind + 1) - 1
              fin = beg + nptsy(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine keveny_wf

subroutine kevenz_wf
  ! ====================================
  ! even kinetic energy block propagation
  ! ====================================
  use globaali, only: nregx, nregy, nregz 
  use globaali, only: nptsx, nptsy, nptsz 
  use globaali, only: wf, hwf,ehwf 
  use globaali, only: lohkox, lohkoy, lohkoz
  use globaali, only: dimx, dimy
  use globaali, only: bot, rgt                               ! ### MPI ###
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
 
  implicit none
  include 'mpif.h'

  integer :: ind, beg, fin, i,j,k
  integer :: Ax, Bx, Ay, By, xend, yend

  Ax = 1 
  Bx = my_dimx
  Ay = 1
  By = my_dimy

  ! third dimension 
  if(bot .NE. MPI_Proc_Null) yend = By - 1 
  if(bot .EQ. MPI_Proc_Null) yend = By
  if(rgt .NE. MPI_Proc_Null) xend = Bx - 1 
  if(rgt .EQ. MPI_Proc_Null) xend = Bx
  do i = Ax, Bx
     do j = Ay, By
        beg = nptsz(1)
        fin = beg + nptsz(2) - 1   
        do ind = 2, nregz, 2
           if(i.LE.xend .AND. j.LE.yend) hwf(i,j,beg:fin) = hwf(i,j,beg:fin) + matmul(lohkoz(ind)%ke, wf(i,j,beg:fin) ) 
           ehwf(i,j,beg:fin) = ehwf(i,j,beg:fin) + matmul(lohkoz(ind)%ke, wf(i,j,beg:fin) ) 
           if(ind .LT. nregz - 1) then
              beg = fin + nptsz(ind + 1) - 1
              fin = beg + nptsz(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine kevenz_wf




subroutine xdy_wf
  ! ====================================
  ! xdy operate on wf
  ! ====================================
  use globaali, only: nregy, nptsy, wf, lohkox, lohkoy, pwx, hwf, ehwf
  use globaali, only: dimx, dimz
  use globaali, only: rgt                                    ! ### MPI ###
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  include 'mpif.h'

  integer :: ind, beg, fin, i, k, myi
  integer :: Ax, Bx, Ay, xend

  Ax = 1 
  Bx = my_dimx
  Ay = 1


  ! along x dimension 
  if(rgt .NE. MPI_Proc_Null) xend = Bx - 1 
  if(rgt .EQ. MPI_Proc_Null) xend = Bx
  myi = x_begind(my_xind) - 1
  do i = Ax, Bx
     myi = myi + 1 
     ! odd x*dy block
     do k = 1, dimz
        beg =                     Ay 
        fin = beg        + nptsy( y_begblk(my_yind) ) - 1  
        do ind = y_begblk(my_yind), y_endblk(my_yind) - 1, 2
         if(i.LE.xend) hwf(i,beg:fin,k) = hwf(i,beg:fin,k) + matmul(pwx(myi,1) * lohkoy(ind)%dy, wf(i,beg:fin,k) ) 
           ehwf(i,beg:fin,k) =  ehwf(i,beg:fin,k) + matmul(pwx(myi,1) * lohkoy(ind)%dy, wf(i,beg:fin,k) ) 
           if(ind .LT. nregy - 1) then
              beg = fin + nptsy(ind + 1) - 1
              fin = beg + nptsy(ind + 2) - 1
           endif
        enddo
     enddo
     ! even x*dy block
     do k = 1, dimz
        beg = Ay + nptsy(y_begblk(my_yind)     ) - 1
        fin = beg               + nptsy(y_begblk(my_yind) + 1 ) - 1   
        do ind = y_begblk(my_yind) + 1, y_endblk(my_yind), 2
           if(i.LE.xend) hwf(i,beg:fin,k) = hwf(i,beg:fin,k) + matmul(pwx(myi,1) * lohkoy(ind)%dy, wf(i,beg:fin,k) )  
           ehwf(i,beg:fin,k) = ehwf(i,beg:fin,k) + matmul(pwx(myi,1) * lohkoy(ind)%dy, wf(i,beg:fin,k) ) 
           if(ind .LT. nregy - 1) then
              beg = fin + nptsy(ind + 1) - 1
              fin = beg + nptsy(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine xdy_wf




subroutine ydx_wf
  ! ====================================
  ! ydx operate on wf
  ! ====================================
  use globaali, only: nregx, nptsx, wf, lohkox, lohkoy, pwy, hwf,ehwf
  use globaali, only: dimy, dimz
  use globaali, only: bot                                    ! ### MPI ###
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  include 'mpif.h'

  integer :: ind, beg, fin, j, k, myj
  integer :: Ax, Ay, By, yend

  Ax = 1 
  Ay = 1
  By = my_dimy


  ! along y dimension 
  if(bot .NE. MPI_Proc_Null) yend = By - 1 
  if(bot .EQ. MPI_Proc_Null) yend = By
  myj = y_begind(my_yind) - 1
  do j = Ay, By
     myj = myj + 1 
     ! odd y*dx block
     do k = 1, dimz
        beg = Ax 
        fin = beg + nptsx( x_begblk(my_xind) ) - 1  
        do ind = x_begblk(my_xind), x_endblk(my_xind) - 1, 2
           if(j.LE.yend) hwf(beg:fin,j,k) = hwf(beg:fin,j,k) + matmul(pwy(myj,1) * lohkox(ind)%dx, wf(beg:fin,j,k) ) 
           ehwf(beg:fin,j,k) =  ehwf(beg:fin,j,k) + matmul(pwy(myj,1) * lohkox(ind)%dx, wf(beg:fin,j,k) ) 
           if(ind .LT. nregx - 1) then
              beg = fin + nptsx(ind + 1) - 1
              fin = beg + nptsx(ind + 2) - 1
           endif
        enddo
     enddo
     ! even y*dx block
     do k = 1, dimz
        beg = Ax + nptsx(x_begblk(my_xind)    ) - 1
        fin = beg               + nptsx(x_begblk(my_xind) + 1) - 1   
        do ind = x_begblk(my_xind) + 1, x_endblk(my_xind), 2
           if(j.LE.yend) hwf(beg:fin,j,k) = hwf(beg:fin,j,k) + matmul(pwy(myj,1) * lohkox(ind)%dx, wf(beg:fin,j,k) )  
           ehwf(beg:fin,j,k) = ehwf(beg:fin,j,k) + matmul(pwy(myj,1) * lohkox(ind)%dx, wf(beg:fin,j,k) ) 
            if(ind .LT. nregx - 1) then
              beg = fin + nptsx(ind + 1) - 1
              fin = beg + nptsx(ind + 2) - 1
           endif
        enddo
     enddo
  enddo

end subroutine ydx_wf


subroutine oddpz_wf
  ! ====================================
  ! odd linear momentum operation
  ! ====================================
  use globaali, only: nregx, nregy, nregz
  use globaali, only: nptsx, nptsy, nptsz
  use globaali, only: wf, hwf 
  use globaali, only: lohkox, lohkoy, lohkoz
  use globaali, only: dimx, dimy 
  use globaali, only: bot, rgt                               ! ### MPI ###
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
  implicit none
  include 'mpif.h'

  integer :: ind, beg, fin, i,j,k
  integer :: Ax, Bx, Ay, By, xend, yend

  Ax = 1 
  Bx = my_dimx
  Ay = 1
  By = my_dimy


! third dimension 
  if(bot .NE. MPI_Proc_Null) yend = By - 1 
  if(bot .EQ. MPI_Proc_Null) yend = By
  if(rgt .NE. MPI_Proc_Null) xend = Bx - 1 
  if(rgt .EQ. MPI_Proc_Null) xend = Bx
  do i = Ax, xend
     do j = Ay, yend
        beg = 1
        fin = nptsz(1)   
        do ind = 1, nregz - 1, 2
           hwf(i,j,beg:fin) = hwf(i,j,beg:fin) + matmul(lohkoz(ind)%dz, wf(i,j,beg:fin) ) 
           if(ind .LT. nregz - 1) then
              beg = fin + nptsz(ind + 1) - 1
              fin = beg + nptsz(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine oddpz_wf

subroutine evenpz_wf
  ! ====================================
  ! even kinetic energy block propagation
  ! ====================================
  use globaali, only: nregx, nregy, nregz 
  use globaali, only: nptsx, nptsy, nptsz 
  use globaali, only: wf, hwf 
  use globaali, only: lohkox, lohkoy, lohkoz
  use globaali, only: dimx, dimy
  use globaali, only: bot, rgt                               ! ### MPI ###
  use globaali, only: x_begind, x_endind, y_begind, y_endind ! ### MPI ###
  use globaali, only: x_begblk, x_endblk, y_begblk, y_endblk ! ### MPI ###
  use globaali, only: my_xind, my_yind                       ! ### MPI ###
  use globaali, only: my_dimx, my_dimy                       ! ### MPI ###
 
  implicit none
  include 'mpif.h'

  integer :: ind, beg, fin, i,j,k
  integer :: Ax, Bx, Ay, By, xend, yend

  Ax = 1 
  Bx = my_dimx
  Ay = 1
  By = my_dimy

  ! third dimension 
  if(bot .NE. MPI_Proc_Null) yend = By - 1 
  if(bot .EQ. MPI_Proc_Null) yend = By
  if(rgt .NE. MPI_Proc_Null) xend = Bx - 1 
  if(rgt .EQ. MPI_Proc_Null) xend = Bx
  do i = Ax, xend
     do j = Ay, yend
        beg = nptsz(1)
        fin = beg + nptsz(2) - 1   
        do ind = 2, nregz, 2
           hwf(i,j,beg:fin) = hwf(i,j,beg:fin) + matmul(lohkoz(ind)%dz, wf(i,j,beg:fin) ) 
           if(ind .LT. nregz - 1) then
              beg = fin + nptsz(ind + 1) - 1
              fin = beg + nptsz(ind + 2) - 1
           endif
        enddo
     enddo
  enddo
end subroutine evenpz_wf
