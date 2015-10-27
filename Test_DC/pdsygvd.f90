!======================================================================
  subroutine PDSYGVD (n, ia, ja, desch, hmat, ib, jb, descb, b, &
                      iz, jz, descz, z, np, nq, evals)
!======================================================================
!
! pdsygvd computes all the eigenvalues, and optionally, the eigenvectors
! of a real generalized SY-definite eigenproblem, of the form
!
!             sub(A) * x = (lambda) * sub(B) * x
!
! sub(A) denoting A(IA:IA+N-1, JA:JA+N-1) is assumed to be SY
!
! sub(B) denoting B(IB:IB+N-1, JB:JB+N-1) is assumed to be symmetric 
! positive definite.
!----------------------------------------------------------------------

    use precisn, only: wp
    use blacs, only: dlen_, ctxt_, nb_, q, nblock, p_error, mb_

    integer, intent(in)     :: n      ! Hamiltonian dimension
    integer, intent(in)     :: ia, ja, ib, jb, iz, jz
    integer, intent(in)     :: desch(dlen_), descb(dlen_), descz(dlen_)
    integer, intent(in)     :: np                    ! local # H rows
    integer, intent(in)     :: nq                    ! local # H cols

    real(wp), intent(inout) :: hmat(:), b(:), z(:)
    real(wp), intent(out)   :: evals(:)              ! eigenvalues
    character(len=1)        :: jobz, uplo
    integer                 :: ibtype, info, iu
    real(wp), allocatable   :: work(:)
    integer, allocatable    :: iwork(:)
    integer                 :: lwork, liwork
    real(wp), save          :: one = 1.0_wp
    character(len=1)        :: trans, range
    integer                 :: anb, ctxt
    integer                 :: liwmin, lwmin, lwopt, mq0, trilwmin
    integer                 :: mycol, myrow, nb, neig, nn, np0, npcol
    integer                 :: nprow, nps, nq0, status
    real(wp)                :: scale
    integer                 :: idum1(5), idum2(5)

    integer                 :: NUMROC

!----------------------------------------------------------------------
! default parameters:
 
    ibtype =  1      ! A * x = lambda * B * x generalized eigenproblem
    range  = 'A'     ! compute all eigenvalues
    jobz   = 'V'     ! compute both eigenvalues and eigenvectors
    uplo   = 'L'     ! use lower triangles of A and B matrices

!----------------------------------------------------------------------
! Get grid parameters:

    ctxt = desch(CTXT_)
    call BLACS_GRIDINFO (ctxt, nprow, npcol, myrow, mycol)

!----------------------------------------------------------------------
! Test the input parameters:

    info = 0
    if (nprow == -1) then
       info = - (900 + ctxt_)
    else if (desch(ctxt_) /= descb(ctxt_)) then
       info = - (1300 + ctxt_)
    else if (desch(ctxt_) /= descz(ctxt_)) then
       info = - (2600 + ctxt_)
    else
       call CHK1MAT (n, 1, n, 1, ia, ja, desch,  4, info)
       call CHK1MAT (n, 1, n, 1, ib, jb, descb,  8, info)
       call CHK1MAT (n, 1, n, 1, iz, jz, descz, 12, info)

       idum1 = (/1, ICHAR('V'), ICHAR('L'), ICHAR('A'), 1/)
       idum2 = (/1, 2, 3, 4, 5/)
       ! check that distributed matrices are consistent across the grid:
       call PCHK2MAT (n, 1, n, 1, ia, ja, desch, 4, n, 1, n, 1, ib, &
            jb, descb, 8, 5, idum1, idum2, info)
       call PCHK1MAT (n, 1, n, 1, iz, jz, descz, 12, 0, idum1, idum2,&
            info)
    end if

    if (info /= 0) then
       call PXERBLA(ctxt, 'PDSYGVD ', -info)
       return
    end if

!----------------------------------------------------------------------
! Form a Cholesky factorization of sub(B).

    call PDPOTRF (uplo, n, b, ib, jb, descb, info)
    call p_error (info, 'pdpotrf: Cholesky error')
       
!----------------------------------------------------------------------
! Transform problem to standard eigenvalue problem.

    nb = desch(nb_)
    np0 = NUMROC(n, nb, 0, 0, nprow)
    nq0 = NUMROC(n, nb, 0, 0, npcol)
    lwork = 2 * np0 * nb + nq0 * nb + nb * nb
    allocate (work(lwork), stat=status)
    call p_error (status, 'pdsyngst: memory allocation')

    call PDSYNGST (ibtype, uplo, n, hmat, ia, ja, desch, b, ib, jb, &
         descb, scale, work, lwork, info)

    call p_error (info, 'pdsyngst: transform error')
    deallocate (work)

! Solve standard eigenvalue problem.

    nblock = desch(mb_)
    liwork = 7*n + 8*q + 2
    trilwmin = 3*n + MAX(nblock*(np+1), 3*nblock)
    lwork = MAX(1+6*n+2*np*nq, trilwmin) + 2*n
    allocate (work(lwork), iwork(liwork), stat=status)
    call p_error (status, 'pdsyevd: memory allocation')

    call BLACS_BARRIER (ctxt, 'all')
    call PDSYEVD (jobz, uplo, n, hmat, ia, ja, desch, evals, z, iz, jz, &
         descz, work, lwork, iwork, liwork, info)
    call p_error (info, 'pdsygvd: pdsyevd error')

    neig = n

! For sub(A)*x=(ambda)*sub(B)*x backtransform eigenvectors: x = inv(L)'*y

    trans = 'T'

! call PBlas level 3 routine:

    one = 1.0_wp
    call PDTRSM ('Left', uplo, trans, 'Non-unit', n, neig, one, &
         b, ib, jb, descb, z, iz, jz, descz)
    if (scale /= 1.0_wp) call DSCAL(n, scale, evals, 1)

  end subroutine PDSYGVD
