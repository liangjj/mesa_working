PROGRAM SINGLEELECTRON
  !     Calculate regular Coulomb wave function
  USE nrtype
  USE pfedvrmod
  USE singleel
  USE dvrmod
  USE f95_lapack
  USE blas95
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: hamiltonian, Vatom, rebuilt_hamiltonian, tmparray
  REAL(DP), DIMENSION(:), ALLOCATABLE :: work
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: ipivot
  INTEGER(I4B) :: i, j, ii, idim, nmax, l, ll1, NR, lwork, ir1, ir2
  INTEGER(I4B) :: i1start, i1end, ireg, info
  !----------Boundary condition at r=0 or not------------------------------
  logical(LGT) :: r0BC = .false.
  !-----------Complete spatial grid and mapping (for "ndim" dimensions)------------------

  REAL(DP) , dimension(:,:), ALLOCATABLE :: grid
  REAL(DP) , dimension(:,:), ALLOCATABLE :: factor 
  INTEGER(I4B), dimension(:), ALLOCATABLE :: regstart, regstop, ngridstart, ngridstop, ngstart, ngstop
  REAL(DP), DIMENSION(:), ALLOCATABLE :: bdstart, bdstop, bdslope, numfunc
  INTEGER(I4B), dimension(:,:,:), ALLOCATABLE :: nindex
  !-------------Total number of grid points for each dimension--------------- 
  integer , dimension(:), ALLOCATABLE :: ntot

  REAL(DP) :: temp, temp1, temp2
  
  !write(0,*) 'huuuhuuu - 00010'

  ALLOCATE(  num_reg(ndim) , ngridstart(ndim), ngridstop(ndim), ngstart(ndim), ngstop(ndim) )	
  ALLOCATE(  bdstart(ndim) , bdstop(ndim), bdslope(ndim), numfunc(ndim) )

  !write(0,*) 'huuuhuuu - 00020'

  open(20,file='singleparams.inp',status='old')
  read(20,*) num_reg(:)
  read(20,*) numfunc(:)
  read(20,*) bdstart(:)
  read(20,*) bdstop(:)
  read(20,*) bdslope(:)
  read(20,*) CoulLMax
  close(20)
  !================================================================
  nmax=maxval(num_reg(:))
  
  !write(0,*) 'huuuhuuu - 00030'
  
  ALLOCATE( ntot(ndim), regstart(ndim), regstop(ndim) )	
  ALLOCATE( num_fun(ndim,nmax) )
  ALLOCATE( bounds(ndim,nmax,1:2) )
  
  ALLOCATE( mat_reg(ndim,nmax) )
  ALLOCATE( nindex(ndim,nmax,500) )
  
  bounds(:,:,:)=0.d0
  num_fun(:,:)=0

  !write(0,*) 'huuuhuuu - 00040'

  !-------------------------------------------------------------
  ! Generating the regions
  !-------------------------------------------------------------

  DO j=1, ndim
     temp=abs(bdstop(j)-bdstart(j))

     i=1
     bounds(j,i,1)=bdstart(j)
     bounds(j,i,2)=bdstart(j)+temp
     num_fun(j,i)=numfunc(j)

     do i=2, num_reg(j)
        bounds(j,i,1)=bounds(j,i-1,2)
        bounds(j,i,2)=bounds(j,i,1)+temp*exp(bdslope(j)*real(i-1,DP)/real(num_reg(j),DP))
        num_fun(j,i)=numfunc(j)
     end do
  END DO

  !write(0,*) 'huuuhuuu - 00050'

  !--------------------------------------------------------
  ! Check if the regions are continuous
  !--------------------------------------------------------

  DO j=1, ndim
     do i=2,num_reg(j)-1
        if(dabs(bounds(j,i,1) - bounds(j,i-1,2)).gt.1.d-10) then
           write(*,*)
           write(*,*)'  The ',j,'th dimension has problem:'
           write(*,*)'  Gap between element',i-1,'and',i,' program halted'
           write(*,*)' bounds(j,i-1,2)=',bounds(j,i-1,2)
           write(*,*)' bounds(j,i,1)=',bounds(j,i,1)
           write(*,*)
           stop
        end if
     end do
  END DO
  
  !write(0,*) 'huuuhuuu - 00060'

  !--------------------------------------------------------
  ! Find total number of grid points
  ! this is the sum of number of functions in each region
  ! minus the number of overlap points
  !--------------------------------------------------------
  do j=1, ndim
     ntot(j)=sum(num_fun(j,:))-(num_reg(j)-1)
     
     !-----------throw out the first point----------------
     !     ntot(j)=ntot(j)-1
  end do

  NR = maxval(ntot(:))

  !write(0,*) 'huuuhuuu - 00070'

  !-------------------------------------------------------------
  !-------------------------------------------------------------

  ALLOCATE ( factor(ndim,NR) )

  !-------------------------------------------------------------
  !   Allocating matrix array for calculations in each element!
  !-------------------------------------------------------------

  !*********first checking if the number of basis functions are the same for each dimension******
  do i=2, num_reg(1)
     if( num_fun(1,i).ne.num_fun(1,i-1) )then
        write(*,*)'i=',i
        write(*,*)'num_fun(1,i)=',num_fun(1,i),' num_fun(1,i-1)=',num_fun(1,i-1)
        write(*,*)'!!!! WE NEED KEEP THE NUMBER OF BASIS SAME FOR EACH ELEMENT IN THE SAME DIMENSION!!!!'
        STOP
     end if
  end do

  !write(0,*) 'huuuhuuu - 00080'

  !-------------------------------------------------------------
  !-----------throw out the first point for both r1 and r2------
  !-------------------------------------------------------------
  !***********This adjustment is very important for the energy calculation***********
  do i=1, ndim
     ngridstart(i)=1
     ngstart(i)=2
     ngridstop(i)=ntot(i)
     ngstop(i)=ngridstop(i)-1
  end do

  !write(0,*) 'huuuhuuu - 00090'

  !  allocate matrices for kinetic energy
  
  DO j=1, ndim
     do i=1, num_reg(j)
        ALLOCATE(mat_reg(j,i)%ke_mat(num_fun(j,i),num_fun(j,i)), &
             mat_reg(j,i)%eigvec_mat(num_fun(j,i),num_fun(j,i)), &
             mat_reg(j,i)%pt(num_fun(j,i)),   &
             mat_reg(j,i)%wt(num_fun(j,i)),   &
             mat_reg(j,i)%fac1(num_fun(j,i)),   &
             mat_reg(j,i)%fac2(num_fun(j,i)),   &
             mat_reg(j,i)%df(num_fun(j,i),num_fun(j,i)),   &
             mat_reg(j,i)%ddf(num_fun(j,i),num_fun(j,i)),   &
             mat_reg(j,i)%eigval_mat(num_fun(j,i)))
     end do
  END DO

  !write(0,*) 'huuuhuuu - 00100'

  !--------------------------------------------------------------------------------
  ! Set up finite element dvr to get points, weights, and Kinetic-energy matrice
  !--------------------------------------------------------------------------------
  CALL dvr_setup

  !write(0,*) 'huuuhuuu - 00110'

  ALLOCATE( grid(ndim,NR) )

  CALL setup_grid_mapping(nindex,grid,factor,r0BC)

  !write(0,*) 'huuuhuuu - 00120'

  ALLOCATE( Vatom(2:NR,0:CoulLMax), hamiltonian(2:NR,2:NR), rebuilt_hamiltonian(2:NR,2:NR), tmparray(2:NR,2:NR) )
  ALLOCATE( singleevals(2:NR,0:CoulLMax), singleevecs(2:NR,2:NR,0:CoulLMax) )

  idim = maxloc(ntot(:),1)
  do l = 0, CoulLmax
     ll1 = l*(l+1)
     do ii = 2, NR
        Vatom(ii,l) = - Z/grid(idim,ii) + ll1/(2.d0*grid(idim,ii)**2)
     end do
  end do
  !write(0,*) 'huuuhuuu - 00130'
  !write(0,*) 'NR =', NR

  do ii = 1, NR
     write(7999,'9999(G13.7,2x)') grid(idim,ii)
  end do

  lwork = 256*(NR-1)
  ALLOCATE( ipivot(2:NR), work(256*(lwork)) )

  do l = 0, CoulLmax
     hamiltonian(:,:) = 0.d0

     i1start = 2
     ireg = 1
     i1end = i1start + num_fun(idim,ireg) - 2 ! not -1 in first region because we don't take the starting point
     hamiltonian(i1start:i1end,i1start:i1end) = hamiltonian(i1start:i1end,i1start:i1end) + mat_reg(idim,ireg)%ke_mat(2:,2:)
     i1start = i1end
     do ireg = 2, num_reg(idim)
        i1end = i1start + num_fun(idim,ireg) - 1
        hamiltonian(i1start:i1end,i1start:i1end) = hamiltonian(i1start:i1end,i1start:i1end) + mat_reg(idim,ireg)%ke_mat(:,:)
        i1start = i1end
     end do
     !write(0,*) 'i1end =', i1end
     !write(0,*) 'huuuhuuu - 00140'

     do ii = 2, NR
        hamiltonian(ii,ii) = hamiltonian(ii,ii) + Vatom(ii,l)
     end do

     !remove first grid point - boundary condition!
     !hamiltonian(1,:) = 0.d0
     !hamiltonian(:,1) = 0.d0
     
     !write(0,*) 'huuuhuuu - 00150'

     !write(77,'(A,I3)') 'hamiltonian for l =', l
     !do ii = 1, NR
     !   write(77,'(9999(G13.7,2x))') hamiltonian(:,ii)
     !end do
     !write(77,*) ''; write(77,*) ''

     singleevecs(:,:,l) = hamiltonian(:,:)
     
     CALL LA_SYEVD(singleevecs(:,:,l), singleevals(:,l), JOBZ='V')

     ! remove highest eigenvalue and see how hamiltonian changes
     
     !singleevals(NR-1,l) = 0.d0
     
     !rebuilt_hamiltonian(:,:) = 0.d0
     !do ii = 2, NR
     !   rebuilt_hamiltonian(ii,ii) = singleevals(ii-1,l)
     !end do

     !CALL GEMM(singleevecs(:,:,l), rebuilt_hamiltonian(:,:), tmparray)
     !CALL GEMM(tmparray, singleevecs(:,:,l), rebuilt_hamiltonian(:,:), transb=blas_trans)
     
     do ii = 2, NR
        write(2000+l,'(9999(G10.4,2x))') hamiltonian(ii,:)
        !write(4000+l,'(9999(G10.4,2x))') rebuilt_hamiltonian(ii,:) - hamiltonian(ii,:)
     end do
     
     !write(6,'(A,I3,A,G20.14))') "sum of absolute errors in rebuilt_hamiltonian for l=",l," : ", sum(abs(rebuilt_hamiltonian(:,:)-hamiltonian(:,:)))
     !write(6,'(A,I3,A,G20.14))') "RMS of absolute errors in rebuilt_hamiltonian for l=",l," : ", sqrt(sum((rebuilt_hamiltonian(:,:)-hamiltonian(:,:))**2))/NR
     write(6,'(A,I3,A,G20.14))') "sum of non-symmetry of hamiltonian            for l=",l," : ", sum(abs(transpose(hamiltonian(:,:))-hamiltonian(:,:)))

     !CALL DSYTRF('U',NR-1,hamiltonian,NR-1,ipivot,work,lwork,info)
     !if (info /= 0) then
     !   write(0,*) 'info /= 0 in DSYTRF! stopping.'
     !   write(0,*) 'info = ', info
     !end if
     !CALL DSYTRI('U',NR-1,hamiltonian,NR-1,ipivot,work,info)
     !if (info /= 0) then
     !   write(0,*) 'info /= 0 in DSYTRI! stopping.'
     !   write(0,*) 'info = ', info
     !end if
     ! rebuild lower part of matrix
     !do ir1 = 2, NR
     !   do ir2 = ir1+1, NR
     !      hamiltonian(ir2,ir1) = hamiltonian(ir1,ir2)
     !   end do
     !end do
     !CALL GEMM(hamiltonian,rebuilt_hamiltonian,singleevecs(:,:,l))
     !do ir1 = 1, NR-1
     !   singleevecs(ir1,ir1,l) = singleevecs(ir1,ir1,l) - 1.d0
     !end do
     
     !write(6,'(A,I3,A,G20.14))') "sum of non-symmetry of inverse hamiltonian    for l=",l," : ", sum(abs(transpose(hamiltonian(:,:))-hamiltonian(:,:)))
     !write(6,'(A,I3,A,G20.14))') "sum of (H^-1 H - 1) (should be zero)          for l=",l," : ", sum(abs(singleevecs(:,:,l)))


     !where(abs(rebuilt_hamiltonian) < 1.d-13) 
     !   rebuilt_hamiltonian = 0.d0
     !end where
     !do ii = 2, NR
     !   write(3000+l,'(9999(G10.4,2x))') rebuilt_hamiltonian(ii,:)
     !end do

     !write a first line for plotting purposes
     write(8000+l,'(9999(G20.14,2x))') 0.d0, singleevecs(ii,:,l)*0.d0
     do ii = 2, NR
        write(7000+l,'(9999(G20.14,2x))') singleevals(ii,l), -2.d0/(ii-1+l)**2,&
             &                 (singleevals(ii,l) * (ii-1+l)**2 / 2.d0) + 1.d0
        write(8000+l,'(9999(G20.14,2x))') grid(idim,ii), singleevecs(ii,:,l)*factor(idim,ii)
     end do
  end do
END PROGRAM SINGLEELECTRON
