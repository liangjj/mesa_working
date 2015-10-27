subroutine kanta
  !============================================================
  ! allocate memory for the basis and compute the
  ! quadrature points and weights and
  ! DVR functions and their first and second derivatives
  !=============================================================
  use globaali , only : lohkox, lohkoy, lohkoz
  use globaali , only : nptsx, nptsy, nptsz
  use globaali , only : nregx, nregy, nregz
  use globaali , only : bridgepts, boundariesx, boundariesy, boundariesz
  implicit none
  
  integer :: ind,jnd
  
  print*, "Creating basis..."


  ! ===========================
  ! first degree of freedom
  ! ===========================
  allocate(lohkox(nregx))
  do ind = 1, nregx
     ! allocate memory
     allocate(lohkox(ind)%dummy(nptsx(ind)))
     allocate(lohkox(ind)%pt(nptsx(ind)))
     allocate(lohkox(ind)%wt(nptsx(ind)))
     allocate(lohkox(ind)%p(nptsx(ind),nptsx(ind)))
     allocate(lohkox(ind)%dp(nptsx(ind),nptsx(ind)))
     allocate(lohkox(ind)%ddp(nptsx(ind),nptsx(ind)))
     lohkox(ind)%dummy = 0d0
     lohkox(ind)%pt    = 0d0
     lohkox(ind)%wt    = 0d0
     lohkox(ind)%p     = 0d0
     lohkox(ind)%dp    = 0d0
     lohkox(ind)%ddp   = 0d0
          

     ! compute quadrature points
     call gaussq('legendre', nptsx(ind), 0d0, 0d0, 2, bridgepts, lohkox(ind)%dummy, lohkox(ind)%pt, lohkox(ind)%wt)
     lohkox(ind)%pt  = (boundariesx(ind+1)-boundariesx(ind)) *  (lohkox(ind)%pt + 1) / 2d0 + boundariesx(ind)
     lohkox(ind)%wt  = lohkox(ind)%wt * (boundariesx(ind+1)-boundariesx(ind)) /2d0  ! normalize
     
     ! compute first and second derivatives
     call lgngr(lohkox(ind)%p,lohkox(ind)%dp,lohkox(ind)%ddp,lohkox(ind)%pt,lohkox(ind)%pt,nptsx(ind),nptsx(ind),'')
     
     deallocate(lohkox(ind)%dummy)

  end do



  ! ===========================
  ! second degree of freedom
  ! ===========================
  allocate(lohkoy(nregy))
  do ind = 1, nregy
     ! allocate memory
     allocate(lohkoy(ind)%dummy(nptsy(ind)))
     allocate(lohkoy(ind)%pt(nptsy(ind)))
     allocate(lohkoy(ind)%wt(nptsy(ind)))
     allocate(lohkoy(ind)%p(nptsy(ind),nptsy(ind)))
     allocate(lohkoy(ind)%dp(nptsy(ind),nptsy(ind)))
     allocate(lohkoy(ind)%ddp(nptsy(ind),nptsy(ind)))
     
     ! compute quadrature points
     call gaussq('legendre', nptsy(ind), 0d0, 0d0, 2, bridgepts, lohkoy(ind)%dummy, lohkoy(ind)%pt, lohkoy(ind)%wt)
     lohkoy(ind)%pt  = (boundariesy(ind+1)-boundariesy(ind)) *  (lohkoy(ind)%pt + 1) / 2d0 + boundariesy(ind)
     lohkoy(ind)%wt  = lohkoy(ind)%wt * (boundariesy(ind+1)-boundariesy(ind)) /2d0  ! normalize
     
     ! compute first and second derivatives
     call lgngr(lohkoy(ind)%p,lohkoy(ind)%dp,lohkoy(ind)%ddp,lohkoy(ind)%pt,lohkoy(ind)%pt,nptsy(ind),nptsy(ind),'')
     
     deallocate(lohkoy(ind)%dummy)
  end do



  ! ===========================
  ! third degree of freedom
  ! ===========================
  allocate(lohkoz(nregz))
  do ind = 1, nregz
     ! allocate memory
     allocate(lohkoz(ind)%dummy(nptsz(ind)))
     allocate(lohkoz(ind)%pt(nptsz(ind)))
     allocate(lohkoz(ind)%wt(nptsz(ind)))
     allocate(lohkoz(ind)%p(nptsz(ind),nptsz(ind)))
     allocate(lohkoz(ind)%dp(nptsz(ind),nptsz(ind)))
     allocate(lohkoz(ind)%ddp(nptsz(ind),nptsz(ind)))
     
     ! compute quadrature points
     call gaussq('legendre', nptsz(ind), 0d0, 0d0, 2, bridgepts, lohkoz(ind)%dummy, lohkoz(ind)%pt, lohkoz(ind)%wt)
     lohkoz(ind)%pt  = (boundariesz(ind+1)-boundariesz(ind)) *  (lohkoz(ind)%pt + 1) / 2d0 + boundariesz(ind)
     lohkoz(ind)%wt  = lohkoz(ind)%wt * (boundariesz(ind+1)-boundariesz(ind)) /2d0  ! normalize
     
     ! compute first and second derivatives
     call lgngr(lohkoz(ind)%p,lohkoz(ind)%dp,lohkoz(ind)%ddp,lohkoz(ind)%pt,lohkoz(ind)%pt,nptsz(ind),nptsz(ind),'')
     
     deallocate(lohkoz(ind)%dummy)
  end do
  
end subroutine kanta
