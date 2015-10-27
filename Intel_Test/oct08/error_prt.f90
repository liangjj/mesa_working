module error_prt
! routine for printing allocation/deallocation errors
! Time-stamp: "02/01/25 07:50:56 cjn"
  use io_units, only:fo
  implicit none

  private
  public alloc_error, io_error

contains

  subroutine alloc_error (status, routine, ad)
    integer, intent(in)          :: status     ! error code
    character(len=*), intent(in) :: routine    ! subroutine name
    character(len=1), intent(in) :: ad         ! alloc / dealloc
    if (ad == 'a' .or. ad == 'A') then
       write (fo,'(a,a,i6)') TRIM(routine), ': allocation error, &
            &status = ', status
    else
       write (fo,'(a,a,i6)') TRIM(routine), ': deallocation error, &
            &status = ', status
    end if
    stop
  end subroutine alloc_error

  subroutine io_error (ios, routine, unt, rw)
! print file i/o error
    integer, intent(in)          :: ios        ! error code
    character(len=*), intent(in) :: routine    ! subroutine name
    character(len=1), intent(in) :: rw         ! read / write
    integer, intent(in)          :: unt        ! file unit number
    if (rw == 'r' .or. rw == 'R') then
       write (fo,'(a,i6,a,i4)') TRIM(routine) // ': unit =', unt,  &
            ' read error, iostat = ', ios
    else
       write (fo,'(a,i6,a,i4)') TRIM(routine) // ': unit =', unt,  &
            ' write error, iostat = ', ios
    end if
    stop
  end subroutine io_error

end module error_prt
