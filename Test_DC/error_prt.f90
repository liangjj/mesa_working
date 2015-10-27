!----------------------------------------------------------------------
    Subroutine alloc_error (status, routine, ad)
!----------------------------------------------------------------------
!   print allocation error
!----------------------------------------------------------------------

    use io_units, only:fo

    integer, intent(in)      :: status     ! error code
    character(*), intent(in) :: routine    ! subroutine name
    character(1), intent(in) :: ad         ! alloc / dealloc

    if (ad == 'a' .or. ad == 'A') then
       write (fo,'(a,a,i6)') TRIM(routine), ': allocation error, &
                             &status = ', status
    else
       write (fo,'(a,a,i6)') TRIM(routine), ': deallocation error, &
                             &status = ', status
    end if
    stop
    End subroutine alloc_error

!----------------------------------------------------------------------
    Subroutine io_error (ios, routine, unt, rw)
!----------------------------------------------------------------------
!   print file i/o error
!----------------------------------------------------------------------

    use io_units, only:fo

    integer, intent(in)      :: ios        ! error code
    character(*), intent(in) :: routine    ! subroutine name
    character(1), intent(in) :: rw         ! read / write
    integer, intent(in)      :: unt        ! file unit number

    if (rw == 'r' .or. rw == 'R') then
       write (fo,'(a,i6,a,i4)') TRIM(routine) // ': unit =', unt,  &
            ' read error, iostat = ', ios
    else
       write (fo,'(a,i6,a,i4)') TRIM(routine) // ': unit =', unt,  &
            ' write error, iostat = ', ios
    end if
    stop
    End subroutine io_error


