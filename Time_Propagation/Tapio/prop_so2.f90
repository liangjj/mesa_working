subroutine prop_so2
  use globaali, only : iu, gamma, omega, method, pq, io, outdir,outfile
  use globaali, only : aika, dt, tmax 
  use globaali, only : pwx, pwy, pwz, wf, myy, ang
  use globaali, only : maxiter, tol, err, errold, myy
  use globaali, only : rank, root, cart_comm,ierr
  
  implicit none

  logical         :: fileflag
  integer         :: ind, i, j, k, maxind, saveind
  real*8          :: p, q, muold, aikaold, de
  
  p    = pq(1)
  q    = pq(2)
  ind  = 1
  aika = dt
  err  = 1
  maxind = floor (tmax / abs(dt))
  saveind = 0
  errold  = 1
  aikaold = -abs(dt)
  
  ! ===================================================
  ! start second order propagation
  ! ==================================================

  ! do while(err .GT. tol .AND. ind .LT. maxiter)
  do while(abs(aika) .LT. tmax )

     pq(3) = 1
     ! ======= apply propagators =======     
     call p_vdiag
     if(omega .NE. 0) then
        call p_rotate
        call p_kinetic_p
        call p_rotate_cc
     else
        call p_kinetic_p
     endif
     aika  = aika + dt
     call p_vdiag


     ! =================================================
     ind = ind + 1
     if(imag(dt) .NE. 0 .OR. gamma .NE. 0) call normalize
     
     ! save the full wavefunction every n=maxind:th step
     if(mod(ind, 200) .EQ. 0 ) then
        if(rank .EQ. root) then          
           print*, 'Saving wfn...t=',abs(aika)
        endif
        if(io) call seiv_wfn(saveind)
        if(rank .EQ. root) then          
           print*, 'Saved.'
        endif
     endif
     ! save a slice of wavefunction every n=20:th step
     if(mod(ind, 20) .EQ. 0 ) then
        de = abs(errold - err) / errold / abs(abs(aika)-aikaold)
        saveind = saveind + 1
        muold = myy
        errold = err
        aikaold = abs(aika)
        if(rank .EQ. root) then          
           print*, 'Saving slice...t=',abs(aika)
        endif
        if(io) call seiv_slice(saveind)
        if(rank .EQ. root) then          
           print*, 'Saved.'
        endif 
        call virhenormi
        ! ========== Adaptive time stepping for imaginary time propagation ==============
         if(de .LT. tol .AND. real(dt).EQ.0 .AND. abs(dt) .GT. 1e-6) then
           dt = dt / 2d0
           if(rank .EQ. root) print*, 'NEW TIME STEP: ', abs(dt) 
           if(rank .EQ. root) print*, "Constructing kinetic energy block propagators..."
           call ke_props
           if(rank .EQ. root) print*, "Constructing momentum operator block propagators..."
          call mom_props
        endif
        
        ! root spits to file
        if(rank.EQ. root) then
           inquire(22,EXIST = fileflag)
           if(fileflag) then
              open(unit=22,file= outdir // 'DE_' // outfile // '.txt', position='append')
           endif
102        format (E12.5,1X,E12.5)
           write(22,102) de,myy
           close(22)
        endif
        
        ! root spits out
        if(rank .EQ. root) then
           print*, 'Step: ', ind
           print 111, myy, real(ang), err
111        format('Mu: ',F12.5,'  L_z: ',F12.5,'  Err: ',E12.5)
           if(imag(dt) .NE. 0 .AND. myy .GT. 10*muold .AND. abs(aika).GT.11* abs(dt)) then
              print*, 'Propagation diverges. Try smaller time-step.'
              print*, 'Chemical potential: ',    myy
              print*, 'Relative error: ',        err
              print*, 'Number of steps taken: ', ind
              stop
           endif
        endif
     endif
  enddo
  ! root spits out
  if(rank .EQ. root) then
     print*, 'Chemical potential: ',    myy
     print*, 'Relative error: ',        err
     print*, 'Number of steps taken: ', ind
  endif
  !===================================================
  ! end fourth order propagation
  ! ==================================================
     
end subroutine prop_so2




 
