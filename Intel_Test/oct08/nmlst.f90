module nmlst
! Namelist input of essential blacs parameters
! Time-stamp: "2008-10-29 06:35:35 cjn"
  implicit none

  character(LEN=3), save        :: ALSP   ! dataset id

  private
  public rd_nmlst
  public ALSP

contains
  subroutine rd_nmlst
! Read namelist input
    use blacs, only: p, q, ictxt, nblock
    use io_units, only: fo, nfnml
    integer                  :: lsp
    integer                  :: ios, dflag
    namelist /pdgin/ lsp, p, q, nblock, dflag

! default values
    lsp = 1              ! Oleg file suffix
    p = 64;     q = 64   ! BLACS grid dimension
    nblock = 64          ! BLACS block size
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
       call BLACS_ABORT (ictxt, ios)
    end if

    if (dflag > 0) write (fo,pdgin)
  end subroutine rd_nmlst

end module nmlst
