subroutine mom_lohkot
  ! ==========================================
  ! construct momentum operator block matrices 
  ! ==========================================
  
  use globaali, only : io, outdir 
  use globaali, only : nregx, nregy, nregz 
  use globaali, only : nptsx, nptsy, nptsz
  use globaali, only : lohkox, lohkoy,lohkoz
  implicit none
  
  integer :: ind, n, m, Lx, Ly, Lz
  real*8  :: br_elem, sqwn, sqwm
  
  ! =======================================
  ! first dimension 
  ! =======================================  
    do ind = 1, nregx   
     allocate(lohkox(ind)%dx(nptsx(ind),nptsx(ind)))
     lohkox(ind)%dx = 0d0
     Lx   = nptsx(ind)
     do n = 1, Lx 
         do m = 1, Lx 
           ! normalzation factors for function n
           if(n.EQ.1 .AND. ind .GT. 1) then
              sqwm = sqrt(lohkox(ind)%wt(1) + lohkox(ind-1)%wt(nptsx(ind-1)))
           else if(n.EQ.Lx .AND. ind .LT. nregx) then
              sqwm = sqrt(lohkox(ind)%wt(Lx) + lohkox(ind+1)%wt(1))
           else
              sqwm = sqrt(lohkox(ind)%wt(n))
           endif
           ! normalization factors for function m
           if(m.EQ.1 .AND. ind .GT. 1) then
              sqwn = sqrt(lohkox(ind)%wt(1) + lohkox(ind-1)%wt(nptsx(ind-1)))
           else if(m.EQ.Lx .AND. ind .LT. nregx) then
              sqwn = sqrt(lohkox(ind)%wt(Lx) + lohkox(ind+1)%wt(1))
           else
              sqwn = sqrt(lohkox(ind)%wt(m))
           endif
           ! normalized matrix element
           lohkox(ind)%dx(n,m) =  (lohkox(ind)%dx(n,m) + lohkox(ind)%wt(n) * lohkox(ind)%dp(n,m) ) / (sqwn*sqwm)
        enddo
     enddo
  enddo
  
  ! =======================================
  ! second dimension 
  ! =======================================  
  
  do ind = 1, nregy   
     allocate(lohkoy(ind)%dy(nptsy(ind),nptsy(ind)))
     lohkoy(ind)%dy = 0d0
     Ly   = nptsy(ind)
     do n = 1, Ly 
         do m = 1, Ly 
           ! normalzation factors for function n
           if(n.EQ.1 .AND. ind .GT. 1) then
              sqwm = sqrt(lohkoy(ind)%wt(1) + lohkoy(ind-1)%wt(nptsy(ind-1)))
           else if(n.EQ.Ly .AND. ind .LT. nregy) then
              sqwm = sqrt(lohkoy(ind)%wt(Ly) + lohkoy(ind+1)%wt(1))
           else
              sqwm = sqrt(lohkoy(ind)%wt(n))
           endif
           ! normalization factors for function m
           if(m.EQ.1 .AND. ind .GT. 1) then
              sqwn = sqrt(lohkoy(ind)%wt(1) + lohkoy(ind-1)%wt(nptsy(ind-1)))
           else if(m.EQ.Ly .AND. ind .LT. nregy) then
              sqwn = sqrt(lohkoy(ind)%wt(Ly) + lohkoy(ind+1)%wt(1))
           else
              sqwn = sqrt(lohkoy(ind)%wt(m))
           endif
           ! normalized matrix element
           lohkoy(ind)%dy(n,m) = (lohkoy(ind)%dy(n,m) + lohkoy(ind)%wt(n) * lohkoy(ind)%dp(n,m) ) / (sqwn*sqwm)
        enddo
     enddo
  enddo
 
  ! =======================================
  ! third dimension 
  ! =======================================  
  
  do ind = 1, nregz   
     allocate(lohkoz(ind)%dz(nptsz(ind),nptsz(ind)))
     lohkoz(ind)%dz = 0d0
     Lz   = nptsz(ind)
     do n = 1, Lz 
         do m = 1, Lz 
           ! normalzation factors for function n
           if(n.EQ.1 .AND. ind .GT. 1) then
              sqwm = sqrt(lohkoz(ind)%wt(1) + lohkoz(ind-1)%wt(nptsz(ind-1)))
           else if(n.EQ.Lz .AND. ind .LT. nregz) then
              sqwm = sqrt(lohkoz(ind)%wt(Lz) + lohkoz(ind+1)%wt(1))
           else
              sqwm = sqrt(lohkoz(ind)%wt(n))
           endif
           ! normalization factors for function m
           if(m.EQ.1 .AND. ind .GT. 1) then
              sqwn = sqrt(lohkoz(ind)%wt(1) + lohkoz(ind-1)%wt(nptsz(ind-1)))
           else if(m.EQ.Lz .AND. ind .LT. nregz) then
              sqwn = sqrt(lohkoz(ind)%wt(Lz) + lohkoz(ind+1)%wt(1))
           else
              sqwn = sqrt(lohkoz(ind)%wt(m))
           endif
           ! normalized matrix element
           lohkoz(ind)%dz(n,m) = (lohkoz(ind)%dz(n,m) + lohkoz(ind)%wt(n) * lohkoz(ind)%dp(n,m) ) / (sqwn*sqwm)
        enddo
     enddo
  enddo
 

  ! ========================================
  ! normalize bridge elements between blocks
  ! ========================================
  do ind = 1, nregx - 1
     Lx                      = nptsx(ind)
     br_elem                 = ( lohkox(ind)%dx(Lx,Lx) + lohkox(ind+1)%dx(1,1) ) / 2d0
     lohkox(ind)%dx(Lx,Lx)   = br_elem
     lohkox(ind+1)%dx(1,1)   = br_elem 
  enddo
  do ind = 1, nregy - 1
     Ly                      = nptsy(ind)
     br_elem                 = ( lohkoy(ind)%dy(Ly,Ly) + lohkoy(ind+1)%dy(1,1) ) / 2d0
     lohkoy(ind)%dy(Ly,Ly)   = br_elem
     lohkoy(ind+1)%dy(1,1)   = br_elem 
  enddo
  do ind = 1, nregz - 1
     Lz                      = nptsz(ind)
     br_elem                 = ( lohkoz(ind)%dz(Lz,Lz) + lohkoz(ind+1)%dz(1,1) ) / 2d0
     lohkoz(ind)%dz(Lz,Lz)   = br_elem
     lohkoz(ind+1)%dz(1,1)   = br_elem 
  enddo
 
end subroutine mom_lohkot
