subroutine ke_props
  ! ==========================================
  ! construct kinetic energy block propagators
  ! ==========================================
  use globaali, only : gamma, method, pq, iu, dt, io 
  use globaali, only : nregx,nregy,nregz
  use globaali, only : nptsx,nptsy, nptsz
  use globaali, only : lohkox, lohkoy, lohkoz

  implicit none
  
  integer                                  :: j,ind,info
  real*8,     allocatable, dimension(:)    :: work
  real*8,     allocatable, dimension(:)    :: eval
  real*8,     allocatable, dimension(:,:)  :: evec
  complex*16, allocatable, dimension(:,:)  :: matmp
  

  ! ==========================================
  ! first dimension
  ! ==========================================
  do ind = 1, nregx
     allocate(work(3*nptsx(ind)))
     allocate(eval(nptsx(ind)))
     allocate(evec(nptsx(ind),nptsx(ind)))
     allocate(matmp(nptsx(ind),nptsx(ind)))
     if(allocated(lohkox(ind)%kprop) .EQ. 0) then
        allocate(lohkox(ind)%kprop(nptsx(ind),nptsx(ind)))
        if(method .EQ. 'SO4') then
           allocate(lohkox(ind)%kpropq(nptsx(ind),nptsx(ind)))
        endif
     endif


     ! ====================
     ! diagonalize block
     ! ====================
     evec = lohkox(ind)%ke
     call dsyev('v', 'u', nptsx(ind), evec, nptsx(ind), eval, work, 3*nptsx(ind), info) 
     ! ====================
     ! construct propagator
     ! ====================
     matmp = transpose(evec)
     do j = 1,nptsx(ind)
        matmp(j,:) = exp(1/(iu-gamma)*pq(1)*dt*eval(j) / (1 + mod(ind,2)) ) * matmp(j,:) 
     enddo
     lohkox(ind)%kprop = matmul(evec,matmp)
     
     if(method .EQ. 'SO4') then
        matmp = transpose(evec)
        do j = 1,nptsx(ind)
           matmp(j,:) = exp(1/(iu-gamma)*pq(2)*dt*eval(j) / (1 + mod(ind,2)) ) * matmp(j,:) 
        enddo
        lohkox(ind)%kpropq = matmul(evec,matmp)
     endif
     
      
     ! ====================================
     ! check convergence of diagonalization
     ! ====================================
     if(info .NE. 0) then
        print*, 'Info = ',info
     endif
     
     deallocate(matmp)
     deallocate(work)
     deallocate(eval)
     deallocate(evec)
  enddo
  
  ! ==========================================
  ! second dimension
  ! ==========================================
  do ind = 1, nregy
     allocate(work(3*nptsy(ind)))
     allocate(eval(nptsy(ind)))
     allocate(evec(nptsy(ind),nptsy(ind)))
     allocate(matmp(nptsy(ind),nptsy(ind)))
     if(allocated(lohkoy(ind)%kprop) .EQ. 0) then
        allocate(lohkoy(ind)%kprop(nptsy(ind),nptsy(ind)))
        if(method .EQ. 'SO4') then
           allocate(lohkoy(ind)%kpropq(nptsy(ind),nptsy(ind)))
        endif
     endif
     ! ==================
     ! diagonalize block
     ! ==================
     evec = lohkoy(ind)%ke
     call dsyev('v', 'u', nptsy(ind), evec, nptsy(ind), eval, work, 3*nptsy(ind), info) 
     ! ====================
     ! construct propagator
     ! ====================
     matmp = transpose(evec)
     do j = 1,nptsy(ind)
        matmp(j,:) = exp(1/(iu-gamma)*pq(1)*dt*eval(j) / (1 + mod(ind,2)) ) * matmp(j,:) 
     enddo
     lohkoy(ind)%kprop = matmul(evec,matmp)
     
     if(method .EQ. 'SO4') then
        matmp = transpose(evec)
        do j = 1,nptsy(ind)
           matmp(j,:) = exp(1/(iu-gamma)*pq(2)*dt*eval(j) / (1 + mod(ind,2)) ) * matmp(j,:) 
        enddo
        lohkoy(ind)%kpropq = matmul(evec,matmp)
     endif
     
     
     ! ====================================
     ! check convergence of diagonalization
     ! ====================================
     if(info .NE. 0) then
        print*, 'Info = ',info
     endif
     
     deallocate(matmp)
     deallocate(work)
     deallocate(eval)
     deallocate(evec)
  enddo
  
  ! ==========================================
  ! third dimension
  ! ==========================================
  do ind = 1, nregz
     allocate(work(3*nptsz(ind)))
     allocate(eval(nptsz(ind)))
     allocate(evec(nptsz(ind),nptsz(ind)))
     allocate(matmp(nptsz(ind),nptsz(ind)))
     if(allocated(lohkoz(ind)%kprop) .EQ. 0) then
        allocate(lohkoz(ind)%kprop(nptsz(ind),nptsz(ind)))
        if(method .EQ. 'SO4') then
           allocate(lohkoz(ind)%kpropq(nptsz(ind),nptsz(ind)))
        endif
     endif
     
     ! ==================
     ! diagonalize block
     ! ==================
     evec = lohkoz(ind)%ke
     call dsyev('v', 'u', nptsz(ind), evec, nptsz(ind), eval, work, 3*nptsz(ind), info) 
     ! ====================
     ! construct propagator
     ! ====================
     matmp = transpose(evec)
     do j = 1,nptsz(ind)
        matmp(j,:) = exp(1/(iu-gamma)*pq(1)*dt*eval(j) / (1 + mod(ind,2)) ) * matmp(j,:) 
     enddo
     lohkoz(ind)%kprop = matmul(evec,matmp)
     
     if(method .EQ. 'SO4') then
        matmp = transpose(evec)
        do j = 1,nptsz(ind)
           matmp(j,:) = exp(1/(iu-gamma)*pq(2)*dt*eval(j) / (1 + mod(ind,2)) ) * matmp(j,:) 
        enddo
        lohkoz(ind)%kpropq = matmul(evec,matmp)
     endif
     
     
     ! ====================================
     ! check convergence of diagonalization
     ! ====================================
     if(info .NE. 0) then
        print*, 'Info = ',info
     endif
     
     deallocate(matmp)
     deallocate(work)
     deallocate(eval)
     deallocate(evec)
  enddo
  
end subroutine ke_props
