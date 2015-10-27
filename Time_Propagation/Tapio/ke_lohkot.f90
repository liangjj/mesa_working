subroutine ke_lohkot
  ! =======================================
  ! construct kinetic energy block matrices 
  ! =======================================
  use globaali, only : io, outdir 
  use globaali, only : nregx, nregy, nregz 
  use globaali, only : nptsx, nptsy, nptsz 
  use globaali, only : lohkox,lohkoy,lohkoz
  
  implicit none
  
  integer :: ind, n, m, Lx, Ly, Lz
  real*8  :: br_elem, sqwn, sqwm
  
  ! =======================================
  ! first dimension 
  ! =======================================  

  ! to make second derivatives exist at bridge points
  do ind = 1, nregx   
     allocate(lohkox(ind)%ke(nptsx(ind),nptsx(ind)))
     lohkox(ind)%ke = 0d0
     Lx = nptsx(ind)
     do m = 1, Lx 
        ! surface terms ala Bloch
        lohkox(ind)%ke(1,m)  =    lohkox(ind)%dp(1,m)
        lohkox(ind)%ke(Lx,m) =  - lohkox(ind)%dp(Lx,m)
     enddo
  enddo

  do ind = 1, nregx   
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
           lohkox(ind)%ke(n,m) = -0.5d0 * (lohkox(ind)%ke(n,m) + lohkox(ind)%wt(n) * lohkox(ind)%ddp(n,m) ) / (sqwn*sqwm)
        enddo
     enddo
  enddo
  
  ! =======================================
  ! second dimension 
  ! =======================================  
  
  ! to make second derivatives exist at bridge points
  do ind = 1, nregy   
     allocate(lohkoy(ind)%ke(nptsy(ind),nptsy(ind)))
     lohkoy(ind)%ke = 0d0
     Ly = nptsy(ind)
     do m = 1, Ly 
        ! surface terms ala Bloch
        lohkoy(ind)%ke(1,m)  =    lohkoy(ind)%dp(1,m)
        lohkoy(ind)%ke(Ly,m) =  - lohkoy(ind)%dp(Ly,m)
     enddo
  enddo

  do ind = 1, nregy   
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
           lohkoy(ind)%ke(n,m) = -0.5d0 * (lohkoy(ind)%ke(n,m) + lohkoy(ind)%wt(n) * lohkoy(ind)%ddp(n,m) ) / (sqwn*sqwm)
        enddo
     enddo
  enddo
  
  ! =======================================
  ! third dimension 
  ! =======================================  
  
  ! to make second derivatives exist at bridge points
  do ind = 1, nregz   
     allocate(lohkoz(ind)%ke(nptsz(ind),nptsz(ind)))
     lohkoz(ind)%ke = 0d0
     Lz = nptsz(ind)
     do m = 1, Lz 
        ! surface terms ala Bloch
        lohkoz(ind)%ke(1,m)  =    lohkoz(ind)%dp(1,m)
        lohkoz(ind)%ke(Lz,m) =  - lohkoz(ind)%dp(Lz,m)
     enddo
  enddo

  do ind = 1, nregz   
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
           lohkoz(ind)%ke(n,m) = -0.5d0 * (lohkoz(ind)%ke(n,m) + lohkoz(ind)%wt(n) * lohkoz(ind)%ddp(n,m) ) / (sqwn*sqwm)
        enddo
     enddo
  enddo
  

  ! ========================================
  ! normalize bridge elements between blocks
  ! ========================================
  do ind = 1, nregx - 1
     Lx                      = nptsx(ind)
     br_elem                 = ( lohkox(ind)%ke(Lx,Lx) + lohkox(ind+1)%ke(1,1) ) / 2d0
     lohkox(ind)%ke(Lx,Lx)   = br_elem
     lohkox(ind+1)%ke(1,1)   = br_elem 
  enddo
  do ind = 1, nregy - 1
     Ly                      = nptsy(ind)
     br_elem                 = ( lohkoy(ind)%ke(Ly,Ly) + lohkoy(ind+1)%ke(1,1) ) / 2d0
     lohkoy(ind)%ke(Ly,Ly)    = br_elem
     lohkoy(ind+1)%ke(1,1)   = br_elem 
  enddo
  do ind = 1, nregz - 1
     Lz                      = nptsz(ind)
     br_elem                 = ( lohkoz(ind)%ke(Lz,Lz) + lohkoz(ind+1)%ke(1,1) ) / 2d0
     lohkoz(ind)%ke(Lz,Lz)   = br_elem
     lohkoz(ind+1)%ke(1,1)   = br_elem 
  enddo
  
end subroutine ke_lohkot
