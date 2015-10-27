module rd_dcmat
! Time-stamp: "2008-11-02 04:26:21 Clifford Noble"
! Read dc matrix computed by Oleg.
  use io_units, only: fo, nfnml
  use blacs, only: p_error, ctxt, mycol, myrow, dlen_, mb_, ctxt_
  use precisn, only: wp
  implicit none

  real(wp), allocatable, save :: hmat(:) ! distributed Hamiltonian
  real(wp), allocatable, save :: z(:)    ! distributed Eigenvectors
  real(wp), allocatable, save :: b(:)    ! distributed Overlaps
  integer, save               :: desch(dlen_) ! Hamiltonian descriptor
  integer, save               :: descz(dlen_) ! Eigenvector descriptor
  integer, save               :: descb(dlen_) ! Overlap descriptor

  integer, save               :: f = 11  ! file unit number

  private
  public read_dcmat
  public hmat, desch, z, descz, b, descb

contains

  subroutine read_dcmat (desch, AF_int, irread, icread, n, np, nq)
! Read overlap and interaction matrix files produced
! by Oleg and distribute them across the processor grid
! ready for scalapack diagonalization
    use nmlst, only: ALSP
    integer, intent(in)    :: desch(:)       ! interaction descriptor
    character(len=*), intent(in) :: AF_int   ! oleg's file name
    integer, intent(in)    :: irread, icread ! coords of io-processor
    integer, intent(out)   :: n              ! H dimension
    integer, intent(out)   :: np, nq         ! local # H rows, cols
    real(wp), allocatable  :: work(:,:)      ! local work matrix
    logical                :: ioprocessor    ! io node flag
    real(wp)               :: alpha, beta, s
    integer                :: descwork(dlen_)
    integer                :: status, lwork, ios, ctxt, nprow, npcol,&
         myrow, mycol, i, ldh, rsrc, csrc, mm, nn, mb, nb, ldw,      &
         j, ich, jch, idim, jdim, ishift, jshift, n0, n1, n2, nhm,   &
         mhm, istart, jstart, info, ns
    integer                :: numroc
    integer, allocatable   :: ns2(:), ns1(:), iprm(:)
    integer                :: parm(1), parms(6)
    character(len=20)      :: AF

    lwork = desch(mb_)   ! blocking factor for distributed matrices
    ctxt = desch(ctxt_)  ! blacs context

    call BLACS_GRIDINFO (ctxt, nprow, npcol, myrow, mycol)
    ioprocessor = (myrow == irread .and. mycol == icread)
 
     if (ioprocessor) then
       i = INDEX(AF_int,'.'); AF = AF_int(1:i) // ALSP
       open (unit=f, file=TRIM(AF), form='unformatted', status='old',&
            access='sequential', action='read', iostat=ios,    &
            position='rewind')
!      iomsg = msg) 
       if (ios /= 0) call p_error (ios,'read_dcmat: file open error')

! n2 = # double continuum channels
! n1 = # single continuum channels
! n0 = # bound configurations
       read (f) n2, n1, n0
       
! ns2(1:n2) = pointer to separate channels in vector p
! ns1(1:n1) = pointer to separate channels in vector f
! ns = # B-splines
! nhm = dim of interaction matrix (after deleting B-splines)
! mhm = dim of interaction matrix (before deleting B-splines)
       if (allocated(ns2)) then
          deallocate(ns2,stat=status)
          if (status /= 0) call p_error (status,'read_dcmat: &
               &deallocation')
       end if
       if (allocated(ns1)) then
          deallocate(ns1,stat=status)
          if (status /= 0) call p_error (status,'read_dcmat: &
               &deallocation')
       end if
       allocate (ns2(0:n2), ns1(0:n1), stat=status)
       if (status /= 0) call p_error (status,'read_dcmat: ns2 &
            &allocation')

       read (f) ns2, ns1, ns, mhm, nhm
       write (*,*) 'nhm = ', nhm, ' ns = ', ns, ' mhm = ', mhm

! iprm = pointer to deleted B-splines
       if (allocated(iprm)) then
          deallocate (iprm, stat=status)
          if (status /= 0) call p_error (status,'read_dcmat: iprm &
               &deallocation')
       end if
       allocate (iprm(mhm), stat=status)
       if (status /= 0) call p_error (status,'read_dcmat: iprm &
            &allocation')
       read (f) iprm(1:mhm)

! delete unwanted arrays:
       deallocate (iprm, ns2, ns1, stat=status)
       if (status /= 0) call p_error (status, 'read_dcmat: &
            &deallocation')

! broadcast dimension data to all processes:
       parm = (/ nhm /)
       call igebs2d (ctxt, 'ALL', ' ', 1, 1, parm, 1)
    else
       call igebr2d (ctxt, 'ALL', ' ', 1, 1, parm, 1, irread, icread)
       nhm = parm(1)
    end if
    n = nhm    ! pass back to rest of code

! initialize array descriptor for the overlap and interaction matrices:
! np, nq are matrix size on current processor
    rsrc = 0;    csrc = 0
    np = numroc (nhm, lwork, myrow, rsrc, nprow)   ! local row dimension
    nq = numroc (nhm, lwork, mycol, csrc, npcol)   ! local col dimension
    ldh = MAX(1, np)

    call descinit (desch, nhm, nhm, lwork, lwork, rsrc, csrc, ctxt, &
         ldh, info)
    if (info /= 0) call p_error (info, 'read_dcmat: descriptor error')
    descz = desch  ! evectors have same description vector as h
    descb = desch  ! overlap has same descriptor as h
    call BLACS_BARRIER (ctxt, 'all')

! Allocate local storage:
    allocate (hmat(np*nq), b(np*nq), z(np*nq), stat=status)
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

! Now repeat for the
! interaction matrix for two-channels submatrix

    vblocks: do
       if (ioprocessor) then
          read (f,iostat=ios) ich, jch, ishift, jshift, idim, jdim
          if (ios == 0 .or. ios == 219) then
             parms = (/ich, jch, ishift, jshift, idim, jdim/)
             call igebs2d (ctxt, 'ALL', ' ', 6, 1, parms, 6)
          else
             call p_error (ios,'dcmat_read: v read')
          end if
       else
          call igebr2d (ctxt, 'ALL', ' ', 6, 1, parms, 6, irread, &
               icread)
          ich = parms(1);      jch = parms(2);    ishift = parms(3)
          jshift = parms(4);  idim = parms(5);      jdim = parms(6)
       end if
       if (ich <= 0) exit vblocks

       if (allocated(work)) then
          deallocate (work, stat=status)
          if (status /= 0) call p_error (status, 'dist: work-&
               &deallocation')
       end if
       allocate (work(idim,jdim), stat=status)
       if (status /= 0) call p_error (status, 'dist: work-allocation')
       work = 0.0_wp
! define blocking factors for work:
       mm = idim;      nn = jdim
       mb = mm;        nb = nn
! source processor
       rsrc = irread;  csrc = icread
       ldw = idim
       call DESCSET (descwork, mm, nn, mb, nb, rsrc, csrc, ctxt, ldw)

       if (ioprocessor) then
          vblock: do 
             read (f) i, j, S
             if (i <= 0) exit vblock
             work(i,j) = s
          end do vblock
       end if
! add block into distributed H matrix:
       alpha = 1.0_wp;   beta = 0.0_wp
       call pdgeadd ('n', idim, jdim, alpha, work, 1, 1, &
            descwork, beta, hmat, ishift+1, jshift+1, desch)
    end do vblocks
    call BLACS_BARRIER (ctxt, 'all')

! now clean up
    deallocate (work, stat=status)
    if (status /= 0) call p_error (status, 'dist: work-deallocation')
  end subroutine read_dcmat

end module rd_dcmat
