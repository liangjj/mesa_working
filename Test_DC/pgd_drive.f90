!======================================================================
    subroutine pgdiag 
!======================================================================
! control the distribution, diagonalization and output
!----------------------------------------------------------------------

    use precisn, only: wp
    use io_units, only: nuh
    use blacs, only: ctxt, nblock, ctxt_, p_error, nb_, mb_, &
        csrc_, p, q, myrow, mycol, io_processor, p0, q0
    use dc_matrix

    integer                  :: i, status, np, nq, n
    real(wp)                 :: alpha, beta

!----------------------------------------------------------------------
! distribute input matrices:

    call read_dcmat (n, np, nq)

    call BLACS_BARRIER (ctxt, 'all')   ! temporary precaution

!----------------------------------------------------------------------
! allocate space to hold the eigenvalues

    if (allocated(eval)) then
     deallocate (eval, stat=status)
     if (status /= 0) call p_error (status,'pgdiag: eval deallocation')
    end if
    allocate (eval(n), stat=status)
    if (status /= 0) call p_error (status,'pgdiag: eval allocation')

!----------------------------------------------------------------------
! initialize distributed eigenvector array, z.

    alpha = 0.0_wp       ! off-diagonal elements
    beta = 1.0_wp        ! diagonal elements
    call pdlaset ('all', n, n, alpha, beta, z, 1, 1, descz)

    call BLACS_BARRIER (ctxt, 'all')    ! precautionary

!----------------------------------------------------------------------
! parallel solution of generalized eigenproblem....
    
    call PDSYGVD1 (n, 1, 1, desch, hmat, 1, 1, descb, b, &
         1, 1, descz, z, np, nq, eval)
    call BLACS_BARRIER (ctxt, 'all')    ! precautionary

!----------------------------------------------------------------------
! write eigenvectors:

    call write_evec (n)

    deallocate (eval, stat=status)

    if (status /= 0) call p_error (status, 'pgdiag: eval &
         &deallocation')

    call BLACS_BARRIER (ctxt, 'all')

  End subroutine pgdiag

