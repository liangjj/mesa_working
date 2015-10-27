module nmlst
! Namelist input of essential blacs parameters
! Time-stamp: "2008-11-02 04:12:49 Clifford Noble"
  implicit none

  character(LEN=3), save        :: ALSP   ! dataset id

  private
  public rd_nmlst
  public ALSP

contains
  subroutine rd_nmlst
! Read namelist input
    use blacs, only: p, q, ictxt, nblock, p_error
    use io_units, only: fo, nfnml
    integer                  :: lsp
    integer                  :: ios, dflag
    namelist /pdgin/ lsp, p, q, nblock, dflag

! default values
    lsp = 1              ! Oleg file suffix
    p = 8;     q = 8     ! BLACS grid dimension
    nblock = 16          ! BLACS block size
    dflag = 0            ! print flag

    open (unit=nfnml, file = 'pgdnl.dat', status='old', &
         position = 'rewind', form='formatted', iostat= ios)
    if (ios /= 0) call BLACS_ABORT (ictxt, ios)

    read (nfnml, pdgin, iostat=ios)
    if (ios /= 0) then
       write (fo, pdgin)
       call BLACS_ABORT (ictxt, ios)
    end if

    close (nfnml)

    write (ALSP,'(i3.3)') lsp   ! convert to characters

    if (p <= 0 .OR. q <= 0) then ! invalid grid dimensions
       write (fo,pdgin)
       call p_error (p+q, 'nmlst: invalid BLACS grid')
    end if

    if (dflag > 0) write (fo,pdgin)
  end subroutine rd_nmlst

end module nmlst
