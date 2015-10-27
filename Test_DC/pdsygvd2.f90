!======================================================================
  subroutine PDSYGVD2 (n, ia, ja, desch, hmat, ib, jb, descb, b, &
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
    real(wp), allocatable   :: work(:), GAP(:)
    integer, allocatable    :: iwork(:), iclustr(:)
    integer                 :: lwork, liwork
    real(wp), save          :: one = 1.0_wp
    character(len=1)        :: trans, range
    integer                 :: anb, ctxt
    integer                 :: liwmin, lwmin, lwopt, mq0, trilwmin
    integer                 :: mycol, myrow, nb, neig, nn, np0, npcol
    integer                 :: nprow, nps, nq0, status
    real(wp)                :: scale, ABSTOL
    integer                 :: NUMROC

!----------------------------------------------------------------------
! default parameters:
 
    ibtype =  1      ! A * x = lambda * B * x generalized eigenproblem
    range  = 'I'     ! compute all eigenvalues
    jobz   = 'V'     ! compute both eigenvalues and eigenvectors
    uplo   = 'L'     ! use lower triangles of A and B matrices

!----------------------------------------------------------------------
! Get grid parameters:

    ctxt = desch(CTXT_)
    call BLACS_GRIDINFO (ctxt, nprow, npcol, myrow, mycol)

!----------------------------------------------------------------------
! Transform problem to standard eigenvalue problem.

    nb = desch(nb_)
    np0 = NUMROC(n, nb, 0, 0, nprow)
    nq0 = NUMROC(n, nb, 0, 0, npcol)
    lwork = 5*n + max(5*n,np0*n+2*nb*nb) + ICEIL(n,nprow*npcol)*n
    allocate (work(lwork), GAP(nprow*npcol), stat=status)
    call p_error (status, 'pdsyngst: memory allocation')
    liwork = 6*MAX(n,nprow*npcol+1,4)
    allocate (iwork(liwork), iclustr(2*nprow*npcol), stat=status)
    call p_error (status, 'pdsyngst: memory allocation')

    VL = -10.0_wp; VU=1000000.0_wp;
    IL = 1; IU = n

    ABSTOL = PDLAMCH(ctxt,'U')

    call PDSYGVX (ibtype,jobz,range,uplo,n,hmat,1,1,desch,b,1,1,descb, &
        VL,VU, IL,IU, ABSTOL, m, nz, evals, -one, z,1,1,descz, work, lwork,&
        iwork, ifail, iclustr,GAP, info)

    call p_error (info, 'pdsyngst: transform error')
    deallocate (work,iwork,iclustr,GAP)


  end subroutine PDSYGVD2
