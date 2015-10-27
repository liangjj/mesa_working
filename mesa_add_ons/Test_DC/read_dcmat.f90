!======================================================================
   subroutine read_dcmat (nhm, np, nq)
!======================================================================
! Read overlap and interaction matrix files produced
! by Oleg and distribute them across the processor grid
! ready for scalapack diagonalization
!----------------------------------------------------------------------

    use blacs
    use dc_matrix
    use io_units

    integer, intent(out)   :: np, nq         ! local # H rows, cols
    real(wp), allocatable  :: work(:,:)      ! local work matrix
    real(wp)               :: alpha, beta, s
    integer                :: descwork(dlen_)
    integer                :: status, nprow, npcol,&
         i, ldh, rsrc, csrc, mm, nn, mb, nb, ldw,      &
         j, ich, jch, idim, jdim, ishift, jshift, n0, n1, n2, nhm,   &
         mhm, info, ns,nw
    integer                :: numroc
    integer, allocatable   :: ns2(:), ns1(:), iprm(:)
    integer                :: parms(5)

    call BLACS_PINFO (iam, nprocs) ! find process #, total # processors

    call BLACS_GRIDINFO (ctxt, nprow, npcol, myrow, mycol)
    io_processor = (myrow == p0 .and. mycol == q0)

!---------------------------------------------------------------------- 
! ... read and broadcast dimenshions

    if (io_processor) then
       i = INDEX(AF_int,'.'); write(AF_int(i+1:i+3),'(i3.3)') klsp
       call check_file (AF_int)
       open (nuh, file=AF_int, form='unformatted')
       ! n2 = # double continuum channels
       ! n1 = # single continuum channels
       ! n0 = # bound configurations

       read (nuh) n2, n1, n0
       ! ns2(0:n2) = pointer to separate channels in vector p
       ! ns1(0:n1) = pointer to separate channels in vector f
       ! ns = # B-splines
       ! nhm = dim of interaction matrix (after deleting B-splines)
       ! mhm = dim of interaction matrix (before deleting B-splines)

       if (allocated(ns2)) deallocate(ns2); allocate (ns2(0:n2))
       if (allocated(ns1)) deallocate(ns1); allocate (ns1(0:n1))
       read (nuh) ns2, ns1, ns, mhm, nhm
       write(*,*) 'ns,mhm,nhm=',ns,mhm,nhm
       ! iprm = pointer to deleted B-splines

       if (allocated(iprm)) deallocate (iprm);  allocate (iprm(mhm))
       read (nuh) iprm(1:mhm)

       ! broadcast dimension data to all processes:

       parms = (/nhm, mhm, ns,nhm,mhm/)
       call igebs2d (ctxt, 'ALL', ' ', 5, 1, parms, 5)
    else
       call igebr2d (ctxt, 'ALL', ' ', 5, 1, parms, 5, p0, q0)
       nhm = parms(1); mhm = parms(2); ns = parms(3)
    end if
    call BLACS_BARRIER (ctxt, 'all')

!---------------------------------------------------------------------- 
! initialize array descriptor for the overlap and interaction matrices:
! np, nq are matrix size on current processor

    rsrc = 0;  csrc = 0
    np = numroc (nhm, nblock, myrow, rsrc, nprow)
    nq = numroc (nhm, nblock, mycol, csrc, npcol)
    ldh = MAX(1, np)

    call descinit (desch, nhm, nhm, nblock, nblock, rsrc, csrc, ctxt, &
         ldh, info)
    call p_error (info, 'read_oleg: descriptor error')

    descz = desch  ! evectors have same description vector as h
    descb = desch  ! overlap has same descriptor as h

    ! Allocate local storage:

    allocate (hmat(np*nq));     hmat = 0.0_wp
    allocate (b(np*nq));        b = 0.0_wp
    allocate (z(np*nq));        z = 0.0_wp

!---------------------------------------------------------------------- 
! initialize work array for reading of matrix blocks:

    if (allocated(work)) deallocate (work)
    nw=ns*ns; allocate (work(nw,nw)); work = 0.0_wp
    call DESCSET (descwork, nw, nw, nw, nw, rsrc, csrc, ctxt, nw)
    call BLACS_BARRIER (ctxt, 'all')

!-----------------------------------------------------------------------
! read overlap matrix for two-channels submatrix:
! ich, jch => indices of the block block
! ishift, jshift => block offsets

    do
       if (io_processor) then
          read (nuh) ich, jch,  ishift, jshift,idim, jdim
          parms = (/ich,ishift,jshift,idim,jdim/)
          call igebs2d (ctxt, 'ALL', ' ', 5, 1, parms, 5)
       else
          call igebr2d (ctxt, 'ALL', ' ', 5, 1, parms, 5, p0, q0)
          ich = parms(1)
          ishift = parms(2);   jshift = parms(3)
          idim = parms(4); jdim = parms(5)
       end if
       if (ich <= 0) exit  ! end of blocks 

       work = 0.0_wp

       ! read nonzero block elements:
       if (io_processor) then
          do; read(nuh) i,j,S;  if(i <= 0) exit; work(i,j)=s; end do 
       end if

       ! add block into distributed overlap matrix:
       alpha = 1.0_wp;   beta = 0.0_wp

       istart=ishift+1
       jstart=jshift+1
       call pdgeadd ('notrans', idim, jdim, alpha, work, 1, 1, &
            descwork, beta, b, istart, jstart, descb)

    end do ! over blocks

if (io_processor) write(*,*) 'Overlap done'

    call BLACS_BARRIER (ctxt, 'all')
!-----------------------------------------------------------------------
! repeat for interaction matrix:

    do
       if (io_processor) then
          read (nuh) ich, jch, ishift, jshift, idim, jdim 
          parms = (/ich,ishift,jshift,idim,jdim/)
          call igebs2d (ctxt, 'ALL', ' ', 5, 1, parms, 5)
       else
          call igebr2d (ctxt, 'ALL', ' ', 5, 1, parms, 5, p0, q0)
          ich = parms(1)
          ishift = parms(2);   jshift = parms(3)
          idim = parms(4); jdim = parms(5)
       end if
 
       if (ich <= 0) exit  ! end of blocks 

       work = 0.0_wp

       ! read nonzero block elements:
       if (io_processor) then
          do; read(nuh) i,j,S;  if(i <= 0) exit; work(i,j)=s; end do 
       end if

       ! add block into distributed overlap matrix:
       alpha = 1.0_wp;   beta = 0.0_wp
       call pdgeadd ('notrans', idim, jdim, alpha, work, 1, 1, &
            descwork, beta, hmat, ishift+1, jshift+1, descb)

    end do ! over blocks

if (io_processor) write(*,*) 'Interaction done'

    call BLACS_BARRIER (ctxt, 'all')

!-----------------------------------------------------------------------

    deallocate (work, stat=status)

  end subroutine read_dcmat


