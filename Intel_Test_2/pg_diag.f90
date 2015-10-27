module pg_diag
! Controls the Parallel solution of the Generalized Eigenproblem
! Real, symmetric Hamiltonian, Positive definitive Overlap matrix
! Time-stamp: "2008-11-04 13:58:15 cjn"
  use precisn, only: wp
  use io_units, only: nuh
  implicit none

  real(wp), allocatable, save   :: eval(:) ! eigenvalues

  private
  public pgdiag

contains

  subroutine pgdiag (irread, icread)
! control the distribution, diagonalization and output
    use blacs, only: ctxt, nblock, ctxt_, p_error, nb_, mb_, &
         csrc_, p, q, myrow, mycol, io_processor, p0, q0
    use gen_diag, only: pdsygvd
    use rd_dcmat, only: read_dcmat, hmat, desch, b, descb, z, descz
    integer, intent(in)      :: irread, icread ! io-processor coords
    character(LEN=20)        :: AF_int ! file-name prefix
    integer                  :: i, status, np, nq, n
    real(wp)                 :: alpha, beta

! provide blacs data:
    desch(mb_) = nblock;    desch(ctxt_) = ctxt

! distribute input matrices:
    AF_int = 'dc_mat.'
    call read_dcmat (desch, AF_int, irread, icread, n, np, nq)

    call BLACS_BARRIER (ctxt, 'all')   ! temporary precaution

! allocate space to hold the eigenvalues
    if (allocated(eval)) then
       deallocate (eval, stat=status)
       if (status /= 0) call p_error (status, 'pgdiag: eval &
            &deallocation')
    end if
    allocate (eval(n), stat=status)
    if (status /= 0) call p_error (status,'pgdiag: eval allocation')

! initialize distributed eigenvector array, z.
    alpha = 0.0_wp       ! off-diagonal elements
    beta = 1.0_wp        ! diagonal elements
    call pdlaset ('all', n, n, alpha, beta, z, 1, 1, descz)

    call BLACS_BARRIER (ctxt, 'all')    ! precautionary

! parallel solution of generalized eigenproblem....
    call PDSYGVD (n, 1, 1, desch, hmat, 1, 1, descb, b, &
         1, 1, descz, z, np, nq, eval)
    call BLACS_BARRIER (ctxt, 'all')    ! precautionary

! write eigenvectors:
    call write_evec (n, p0, q0)

    deallocate (eval, stat=status)
    if (status /= 0) call p_error (status, 'pgdiag: eval &
         &deallocation')
    call BLACS_BARRIER (ctxt, 'all')
  end subroutine pgdiag

  subroutine write_evec (n, p0, q0)
! obtain a column of distributed Hamiltonian eigenvectors, z
! local pieces are sent to process with coordinates (p0, q0).
! externals: blacs_gridinfo, descset, pdgeadd, pxerbla
! assume x has dimension n
    use blacs, only: mb_, ctxt_, dlen_, p_error
    use rd_dcmat, only: z, descz
    use nmlst, only: ALSP
    integer, intent(in)   :: n      ! Hamiltonian dimensions
    integer, intent(in)   :: p0, q0 ! i/o proc coords
    character(LEN=20)     :: AF_evc ! Output file prefix
    character(LEN=23)     :: AF     ! output file name
    real(wp), allocatable :: x(:)   ! retrieved eigenvector
    logical               :: ioprocessor
    integer               :: nprow, npcol, myrow, mycol
    integer               :: mm, nn, lwork, ctxt, mb, nb, rsrc, csrc
    integer               :: ldx, isize, istart, iend, i, ios
    integer               :: jsize, jstart, jend, ix, status
    real(wp)              :: alpha, beta
    integer               :: descx(dlen_)

    lwork = descz(mb_)   ! x dimension
    ctxt = descz(ctxt_)  ! context

! get blacs information:
    call BLACS_GRIDINFO (ctxt, nprow, npcol, myrow, mycol)
    ioprocessor = (myrow == p0 .and. mycol == q0)

    if (ioprocessor) then
       AF_evc = 'dc_vecs.'
       i = INDEX(AF_evc,'.'); AF = AF_evc(1:i) // ALSP
       open (unit=nuh, file=TRIM(AF), form='unformatted',         &
            status='unknown', access='sequential', action='write',&
            iostat=ios)
       if (ios /= 0) call p_error (ios,'write_evec: file open error')
       write (nuh) n
       write (nuh) eval
    end if

    allocate (x(n), stat=status)
    if (status /= 0) call p_error (status,'write_evect allocation')
! form blacs descriptor for vector x:
    mm = MAX(1, MIN(n, lwork))    ! normally lwork
    nn = 1
    mb = mm;     nb = nn
    rsrc = p0;   csrc = q0
    ldx = MAX(1, mm)
    call descset (descx, mm, nn, mb, nb, rsrc, csrc, ctxt, ldx)

    alpha = 1.0_wp;     beta = 0.0_wp
    jst: do jstart = 1, n
       jend = jstart
       jsize = 1
       ist: do istart = 1, n, mm
          iend = MIN(n, istart+mm-1)
          isize = iend - istart + 1
          call pdgeadd ('notrans', isize, jsize, alpha, z, istart, &
               jstart, descz, beta, x(istart:iend), 1, 1, descx)
       end do ist
       call BLACS_BARRIER (ctxt, 'all')

       if (ioprocessor) then   
          write (nuh) jstart
          write (nuh) x
          if (jstart <= 3 .or. jstart >= n-2) then
             write (*,*) jstart
             write (*,*) x
          end if
       end if
    end do jst
    deallocate (x, stat=status)
    if (status /= 0) call p_error (status,'write_evect deallocation')
    if (ioprocessor) endfile (unit=nuh)
  end subroutine write_evec
end module pg_diag


