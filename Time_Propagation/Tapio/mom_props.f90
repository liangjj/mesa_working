subroutine mom_props
  ! ==========================================
  ! construct kinetic energy block propagators
  ! ==========================================
  use globaali, only : method, pq, iu, dt, io, omega, gamma
  use globaali, only : pwx, pwy 
  use globaali, only : nregx, nregy 
  use globaali, only : nptsx, nptsy 
  use globaali, only : lohkox, lohkoy 
  use globaali, only : lyprop, lxprop

  implicit none
  
  integer                                      :: j, ind, info, ix, iy, koko, status
  complex*16,     allocatable, dimension(:)    :: work
  complex*16,     allocatable, dimension(:)    :: rwork
  real*8    ,     allocatable, dimension(:)    :: eval
  complex*16,     allocatable, dimension(:,:)  :: evec
  complex*16, allocatable, dimension(:,:)      :: matmp
  
  
  ! ==========================================
  ! first dimension
  ! ==========================================
  if(allocated(lyprop) .EQ. 0)  allocate(lyprop(size(pwx(:,1))))
  do ix = 1, size(pwx(:,1))
     if(allocated(lyprop(ix)%blk) .EQ. 0)  allocate(lyprop(ix)%blk(nregy))
     do ind = 1, nregy
        koko = nptsy(ind)
        allocate(work(3*koko))
        allocate(rwork(3*koko-2))
        allocate(eval(koko))
        allocate(evec(koko,koko))
        allocate(matmp(koko,koko))
        if(allocated(lyprop(ix)%blk(ind)%mprop) .EQ. 0) then
           allocate(lyprop(ix)%blk(ind)%mprop(koko,koko))
           if(method .EQ. 'SO4') then
              allocate(lyprop(ix)%blk(ind)%mpropq(koko,koko))
           endif
        endif

        ! ==================
        ! diagonalize block
        ! ==================
        evec = iu * omega * pwx(ix,1) * lohkoy(ind)%dy  
        call zheev('v', 'u', koko, evec, koko, eval, work, 3*koko, rwork, info) 
        ! ====================
        ! construct propagator
        ! ====================
        matmp = transpose(conjg(evec))
        do j = 1,koko
           matmp(j,:) = exp(1/(iu-gamma)*pq(1)*dt*eval(j) / (2 + 0*mod(ind,2)) ) * matmp(j,:) 
        enddo
        lyprop(ix)%blk(ind)%mprop = matmul(evec,matmp)
        
        if(method .EQ. 'SO4') then
           matmp = transpose(conjg(evec))
           do j = 1,koko
              matmp(j,:) = exp(1/(iu-gamma)*pq(2)*dt*eval(j) / (2 + 0*mod(ind,2)) ) * matmp(j,:) 
           enddo
           lyprop(ix)%blk(ind)%mpropq = matmul(evec,matmp)           
        endif
        

        ! ====================================
        ! check convergence of diagonalization
        ! ====================================
        if(info .NE. 0) then
           print*, 'Info = ',info
        endif
        
        deallocate(matmp)
        deallocate(work)
        deallocate(rwork)
        deallocate(eval)
        deallocate(evec)
        
     enddo
  enddo
  
  
  ! ==========================================
  ! second dimension
  ! ==========================================
  if(allocated(lxprop) .EQ. 0)  allocate(lxprop(size(pwy(:,1))))
  do iy = 1, size(pwy(:,1))
     if(allocated(lxprop(iy)%blk) .EQ. 0) allocate(lxprop(iy)%blk(nregx))
     do ind = 1, nregx
        koko = nptsx(ind)
        allocate(work(3*koko))
        allocate(rwork(3*koko -2))
        allocate(eval(koko))
        allocate(evec(koko,koko))
        allocate(matmp(koko,koko))
        if(allocated(lxprop(iy)%blk(ind)%mprop) .EQ. 0) then
           allocate(lxprop(iy)%blk(ind)%mprop(koko,koko))
           if(method .EQ. 'SO4') then
              allocate(lxprop(iy)%blk(ind)%mpropq(koko,koko))
           endif
        endif
        
        ! ==================
        ! diagonalize block
        ! ==================
        evec = - iu * omega * pwy(iy,1) * lohkox(ind)%dx
        call zheev('v', 'u', koko, evec, koko, eval, work, 3*koko, rwork, info) 
        ! ====================
        ! construct propagator
        ! ====================
        matmp = transpose(conjg(evec))
        do j = 1,koko
           matmp(j,:) = exp(1/(iu-gamma)*pq(1)*dt*eval(j) / (2 + 0*mod(ind,2)) ) * matmp(j,:) 
        enddo
        lxprop(iy)%blk(ind)%mprop = matmul(evec,matmp)
        
        if(method .EQ. 'SO4') then
           matmp = transpose(conjg(evec))
           do j = 1,koko
              matmp(j,:) = exp(1/(iu-gamma)*pq(2)*dt*eval(j) / (2 + 0*mod(ind,2)) ) * matmp(j,:) 
           enddo
           lxprop(iy)%blk(ind)%mpropq = matmul(evec,matmp)
        endif
        
        ! ====================================
        ! check convergence of diagonalization
        ! ====================================
        if(info .NE. 0) then
           print*, 'Info = ',info
        endif
  
        deallocate(matmp)
        deallocate(work)        
        deallocate(rwork)
        deallocate(eval)
        deallocate(evec)
     enddo
  enddo
  
end subroutine mom_props
