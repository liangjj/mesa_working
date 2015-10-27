module trans
  use blacs, only: nprocs, iam, ctxt, nprow, npcol, myrow, mycol, p0, &
       q0, dlen_
  use io_units, only: nud
  implicit none

  private
  public read_dlmat

contains

  subroutine read_dlmat (nhm1, nhm2, np, nq)
! Read dipole matrix files produced by BSR in a B-spline basis
    use nmlst, only: AF_dmat, klsp1, klsp2
    integer, intent(out)   :: nhm1, nhm2     ! dipole matrix dimensions
    integer, intent(out)   :: np, nq         ! local # dipole rows, cols
    real(wp), allocatable  :: work(:,:)      ! local work matrix
    real(wp), allocatable  :: ns2(:), ms2(:)      
    real(wp)               :: alpha, beta, s
    integer                :: descwork(dlen_)
    integer                :: status, nprow, npcol, i, ldd, rsrc, csrc,&
         j, ich, jch, idim, jdim, ishift, jshift, n0, n1, n2, m0, m1,  &
         m2, mhm1, mhm2, info, ns,nw 
    integer                :: numroc
    integer                :: parms(5)

    call BLACS_GRIDINFO (ctxt, nprow, npcol, myrow, mycol)
    io_processor = (myrow == p0 .and. mycol == q0)

! read and broadcast dimenshions
    if (io_processor) then
       i = INDEX(AF_dmat,'.')
       write(AF_dmat(i+1:i+3),'(i3.3)') klsp1
       write(AF_dmat(i+5:i+7),'(i3.3)') klsp2

       open (nud, file=TRIM(AF_dmat), form='unformatted', status='old',&
            iostat=ios, position='rewind', action='read')
       if (ios /= 0) call p_error (ios, 'read_dlmat: AF_dmat open')

       read (nud) n1, n2, n0
       read (nud) m1, m2, m0
       read (nud) ns

       write (*,*) 'n2,n1,n0', n2, n1, n0
       write (*,*) 'm2,m1,m0', m2, m1, m0

       read (nud) mhm1, nhm1
       read (nud) mhm2, nhm2

       write(*,*) 'mhm1, nhm1', mhm1, nhm1
       write(*,*) 'mhm2, nhm2', mhm2, nhm2

! nhm = dim of interaction matrix (after deleting B-splines)
! mhm = dim of interaction matrix (before deleting B-splines)

! broadcast dimension data to all processes:
       parms = (/nhm1, ns, nhm2/)
       call igebs2d (ctxt, 'ALL', ' ', 3, 1, parms, 3)
    else
       call igebr2d (ctxt, 'ALL', ' ', 3, 1, parms, 3, p0, q0)
       nhm1 = parms(1);  ns = parms(2);  nhm2 = parms(3)
    end if
    call BLACS_BARRIER (ctxt, 'all')

! initialize array descriptor for the dipole matrices:
! np, nq are the dipole matrix size on current processor
    rsrc = 0;  csrc = 0
    np = numroc (nhm1, nblock, myrow, rsrc, nprow)
    nq = numroc (nhm2, nblock, mycol, csrc, npcol)
    ldd = MAX(1, np)

    call descinit (descd, nhm1, nhm2, nblock, nblock, rsrc, csrc, ctxt,&
         ldd, info)
    if (info /= 0) call p_error (info, 'read_dlmat: descriptor error')

    ! Allocate local storage:
! Allocate local storage:
    allocate (dip(np*nq), b(np*nq), z(np*nq), stat=status)
    if (status /= 0) call p_error (status,'read_dcmat: hmat memory')
    hmat = 0.0_wp;   b = 0.0_wp
    call BLACS_BARRIER (ctxt, 'all')

! overlap matrix for two-channels submatrix:

    oblocks: do
! ich, jch = indices within block
! ishift, jshift = block offsets
       if (ioprocessor) then
          read (f, iostat=ios) ich, jch, ishift, jshift, idim, jdim
          if (ios == 0 .or. ios == 219) then
             parms = (/ich, jch, ishift, jshift, idim, jdim/)
          else
             call p_error (ios,'read_dcmat: read error')
          end if
          call igebs2d (ctxt, 'ALL', ' ', 6, 1, parms, 6)
       else
          call igebr2d (ctxt, 'ALL', ' ', 6, 1, parms, 6, irread, &
               icread)
          ich = parms(1);      jch = parms(2);     ishift = parms(3); 
          jshift = parms(4);   idim = parms(5);    jdim = parms(6)
       end if
       if (ich <= 0) exit oblocks

! work holds each block of data read from Oleg's file.
! this may have to be resized for each block 
       if (allocated(work)) then
          deallocate (work, stat=status)
          if (status /= 0) call p_error (status, 'dist: work-&
               &deallocation')
       end if
       allocate (work(idim,jdim), stat=status)
       if (status /= 0) call p_error (status, 'dist: work-allocation')
       work = 0.0_wp
! define a pblas descriptor vector for work
! define blocking factors for work:
       mm = idim;      nn = jdim     ! overall dimensions
       mb = mm;        nb = nn       ! local dimensions (only on 1 proc)
       rsrc = irread;  csrc = icread  ! source processor containing work
       ldw = idim                    ! leading dimension of work
       call DESCSET (descwork, mm, nn, mb, nb, rsrc, csrc, ctxt, ldw)

! read nonzero block elements:
       if (ioprocessor) then
          oblock: do 
             read (f) i, j, S
             if (i <= 0) exit oblock
             work(i,j) = s
          end do oblock
       end if
! add block into distributed overlap matrix:
       alpha = 1.0_wp;   beta = 0.0_wp
       call pdgeadd ('n', idim, jdim, alpha, work, 1, 1, &
            descwork, beta, b, ishift+1, jshift+1, descb)
    end do oblocks
    call BLACS_BARRIER (ctxt, 'all')

    allocate (dip(np,nq));     dip = 0.0_wp

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
          read (nud) ich, jch,  ishift, jshift, idim, jdim

!if(ich.gt.0) then
! i=1; if(ich.gt.1) i=ns2(ich-1); idim=ns2(ich)-i
! i=1; if(jch.gt.1) i=ms2(jch-1); jdim=ms2(jch)-i
!end if

          parms = (/ich,ishift,jshift,idim,jdim/)
          call igebs2d (ctxt, 'ALL', ' ', 5, 1, parms, 5)

write(*,*) 'ich,jch',ich,jch,ishift,jshift,idim, jdim

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
          do; read(nud) i,j,S;  if(i <= 0) exit; work(i,j)=s; end do 
       end if

       ! add block into distributed overlap matrix:
       alpha = 1.0_wp;   beta = 0.0_wp

       call pdgeadd ('notrans', idim, jdim, alpha, work, 1, 1, &
            descwork, beta, dip, ishift+1, jshift+1, descd)

    end do ! over blocks

if (io_processor) write(*,*) 'DIP-matrix done'

    deallocate (work, stat=status)

    call BLACS_BARRIER (ctxt, 'all')

  end subroutine read_dlmat


