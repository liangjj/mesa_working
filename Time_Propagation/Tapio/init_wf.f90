subroutine init_wf
! ================================================
! - allocate memory for wavefunction and potential
! - initialize wavefunction, and other necessities
! ================================================

  use globaali, only : method, instate
  use globaali, only : C, pq, io, iu, indir, infile, outdir, outfile 
  use globaali, only : pwx, pwy, pwz
  use globaali, only : wf, hwf, hpsi, locerr, vext, vprop, ehwf
  use globaali, only : lambdax,lambday, lambdaz, mutf
  use globaali, only : dimx, dimy, dimz
! ================= MPI ====================
  use globaali, only : rank, root, cart_comm, ierr
  use globaali, only : x_begind, x_endind, y_begind, y_endind                  
  use globaali, only : my_xind, my_yind                                        
  use globaali, only : my_dimx, my_dimy                                       


  implicit none
  
  integer  :: i, j, k, myi, myj
  real*8   :: pi, rr, ir, rtf2, wfre, wfim
  pi = acos(-1d0)
  
  allocate(wf(my_dimx,my_dimy,dimz))
  allocate(hwf(my_dimx,my_dimy,dimz))
  allocate(ehwf(my_dimx,my_dimy,dimz))
  allocate(hpsi(my_dimx,my_dimy,dimz))
  allocate(vext(my_dimx,my_dimy,dimz))
  allocate(vprop(my_dimx,my_dimy,dimz))
  
  
  if(method .EQ. 'SO4') then
     pq(1) = 1d0 / (4d0 - 4d0**(1d0/3d0))
     pq(2) = 1d0 - 4d0*pq(1)
     if(rank .EQ. root) then                                      
        print*, 'Time propagation method: 4th order split operator'
     endif
  else if(method .EQ. 'SO2') then
     pq(1) = 1d0
     print*, 'Time propagation method: 2nd order split operator'
  else
     print*, 'Chosen time propagation method is unsupported'
     stop
  endif
  
 ! initialize wavefunction
 myi = 0
 myj = 0
 do i = x_begind(my_xind), x_endind(my_xind) 
    myi = myi + 1
    myj = 0 
    do j = y_begind(my_yind), y_endind(my_yind)
        myj = myj + 1
        do k = 1, dimz
           ! =======================================================================================
           if(instate .EQ. 'GA') then      ! GAUSSIAN
              wf(myi,myj,k)     = 1 / (pi**(0.25d0))**3 * exp( -((lambdax*pwx(i,1))**2 + (lambday*pwy(j,1))**2 + (lambdaz*pwz(k,1))**2) / 2d0) * sqrt(pwx(i,2)*pwy(j,2)*pwz(k,2))
              ! ====================================================================================
           else if(instate .EQ. 'TF') then ! THOMAS-FERMI
              rtf2 = ((lambdax*pwx(i,1))**2 + (lambday*pwy(j,1))**2 + (lambdaz*pwz(k,1))**2)
              if(rtf2 .LT. 2*mutf) then     
                 wf(myi,myj,k) =  sqrt(abs(mutf - 0.5*rtf2) / C) * sqrt(pwx(i,2)*pwy(j,2)*pwz(k,2))
              endif
              ! ====================================================================================
           else if(instate .EQ. 'VT') then ! THOMAS-FERMI WITH VORTICITY
              rtf2 = ((lambdax*pwx(i,1))**2 + (lambday*pwy(j,1))**2 + (lambdaz*pwz(k,1))**2)
              if(rtf2 .LT. 2*mutf) then
                 wf(myi,myj,k) =  exp(iu*1*atan2((pwy(j,1) + 0),pwx(i,1))) * sqrt(abs(mutf - 0.5*rtf2) / C) * sqrt(pwx(i,2)*pwy(j,2)*pwz(k,2))
              endif
              ! ====================================================================================
           else if(instate .EQ. 'CC') then ! UNIT
              wf(myi,myj,k) =  1 * sqrt(pwx(i,2)*pwy(j,2)*pwz(k,2))
              ! ==================================================================================== 
           else if(instate .EQ. 'RN') then ! RANDOM
              call random_number(rr)
              call random_number(ir)
              wf(myi,myj,k) =  cmplx(rr,ir) *  sqrt(pwx(i,2)*pwy(j,2)*pwz(k,2))
              ! ====================================================================================
           else
              print*, 'Illegal instate '
              stop
           endif
        enddo
     enddo
  enddo

  ! read from file
  if(instate .EQ. 'IN') then
     open(unit=10,file= indir // infile // '.wfn')
     myi = 0
     myj = 0
     do i = 1, dimx
        if(i.GE.x_begind(my_xind) .AND. i.LE. x_endind(my_xind) ) myi = myi + 1
        myj = 0 
        do j = 1, dimy
           if(j.GE.y_begind(my_yind) .AND. j.LE. y_endind(my_yind) ) myj = myj + 1
           do k = 1, dimz
              read(unit=10,fmt=*) wfre, wfim 
              if(i.GE.x_begind(my_xind) .AND. i.LE. x_endind(my_xind) .AND. j.GE.y_begind(my_yind) .AND. j.LE. y_endind(my_yind) ) then
                 wf(myi,myj,k) = cmplx(wfre, wfim)
                 !if(sqrt(pwx(i,1)**2+pwy(j,1)**2+pwz(k,1)**2) .GT. 20 ) then 
                 !   wf(myi,myj,k)= 0d0
                 !endif
              endif
           enddo
        enddo
     enddo
     close(unit=10)
  endif
  
! make sure initial wf is normalized
  call normalize

! write inital state to file  
  if(io) call seiv_wfn(0)
  if(io) call seiv_slice(0)


end subroutine init_wf
