subroutine ptswts
  ! =========================================================
  ! compute global points and weights  
  ! bridge weights = w_i + w_(i+1)
  ! =========================================================

  use globaali, only : io, outdir 
  use globaali, only : nptsx, nptsy, nptsz
  use globaali, only : nregx, nregy, nregz
  use globaali, only : lohkox, lohkoy, lohkoz
  use globaali, only : pwx, pwy,pwz
  use globaali, only : root, rank
  implicit none
  
  integer :: Lx,Ly,Lz,ind,i,j

  ! ==================================
  ! first dimension
  ! ==================================
  Lx = sum(nptsx(:)) + 1 - nregx
  allocate(pwx(Lx,2))
  pwx(:,2) = 0d0
  j = 1
  do ind = 1, nregx
     pwx(j:j+nptsx(ind)-1,2) = pwx(j:j+nptsx(ind)-1,2) + lohkox(ind)%wt
     do i = 1, nptsx(ind) - 1
        pwx(j,1) = lohkox(ind)%pt(i)
        j = j + 1  
     enddo
  enddo
  pwx(Lx,1) = lohkox(nregx)%pt(nptsx(nregx))



  ! ==================================
  ! second dimension
  ! ==================================
  Ly = sum(nptsy(:)) + 1 - nregy
  allocate(pwy(Ly,2))
  pwy(:,2) = 0d0
  j = 1
  do ind = 1, nregy
     pwy(j:j+nptsy(ind)-1,2) = pwy(j:j+nptsy(ind)-1,2) + lohkoy(ind)%wt
     do i = 1, nptsy(ind) - 1
        pwy(j,1) = lohkoy(ind)%pt(i)
        j = j + 1  
     enddo    
  enddo
  pwy(Ly,1) = lohkoy(nregy)%pt(nptsy(nregy))


  ! ==================================
  ! third dimension
  ! ==================================
  Lz = sum(nptsz(:)) + 1 - nregz
  allocate(pwz(Lz,2))
  pwz(:,2) = 0d0
  j = 1
  do ind = 1, nregz
     pwz(j:j+nptsz(ind)-1,2) = pwz(j:j+nptsz(ind)-1,2) + lohkoz(ind)%wt
     do i = 1, nptsz(ind) - 1
        pwz(j,1) = lohkoz(ind)%pt(i)
        j = j + 1  
     enddo    
  enddo
  pwz(Lz,1) = lohkoz(nregz)%pt(nptsz(nregz))


  ! ==================================
  ! print points and weights to file
  ! ==================================
  if (io .AND. rank .EQ. root) then
     open(unit=11,file= outdir // 'pwx.txt')
     do i = 1,Lx
        write(11,100) pwx(i,1),pwx(i,2)
     enddo
     close(unit=11)
  end if
  if (io .AND. rank .EQ. root) then
     open(unit=11,file= outdir // 'pwy.txt')
     do i = 1,Ly
        write(11,100) pwy(i,1),pwy(i,2)
     enddo
     close(unit=11)
  end if
  if (io .AND. rank .EQ. root) then
     open(unit=11,file= outdir // 'pwz.txt')
     do i = 1,Lz
        write(11,100) pwz(i,1),pwz(i,2)
     enddo
     close(unit=11)
  end if
100 format(E12.5,1X,E12.5/)
  
end subroutine ptswts
