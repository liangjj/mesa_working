!======================================================================
  subroutine write_evec (n)
!======================================================================
! obtain a column of distributed Hamiltonian eigenvectors, z
! local pieces are sent to process with coordinates (p0, q0).
! externals: blacs_gridinfo, descset, pdgeadd, pxerbla
! assume x has dimension n
!----------------------------------------------------------------------

    use io_units
    use blacs
    use dc_matrix

    integer, intent(in)   :: n      ! Hamiltonian dimensions

    real(wp), allocatable :: x(:)   ! retrieved eigenvector

    real(wp)              :: alpha, beta
    integer               :: descx(dlen_)
    integer               :: j
    integer               :: len
    integer               :: lenth
    character(LEN=80)     :: itoc
    character(LEN=80)     :: vector_index

  
! get blacs information:

    call BLACS_GRIDINFO (ctxt, nprow, npcol, myrow, mycol)
    io_processor = (myrow == p0 .and. mycol == q0)

! open output file:

    if (io_processor) then
       i = INDEX(AF_sol,'.'); write(AF_sol(i+1:i+3),'(i3.3)') klsp
       open (nus, file=TRIM(AF_sol), form='unformatted',         &
            status='unknown', access='sequential', action='write',&
            iostat=ios)
!       Call Iosys('open eigen_file as new',0,0,0,'eigen_file')
       if (ios /= 0) call p_error (ios,'write_evec: file open error')
       write (nus) n
       write (nus) eval
!       Call IOsys('write integer matrix_size to eigen_file',1,n,0,' ')
!       Call IOsys('write real eigenvalues to eigen_file',n,eval,0,' ')
       write (fo,*) 'n =',n
       write (fo,'(i10,f16.8)') (i,eval(i),i=1,n)
!       Call IOsys('create real eigenvectors on eigen_file',n*n,0,0,' ')
    end if
    allocate (x(n))

! form blacs descriptor for vector x:

    call descset (descx, n, 1, n, 1, p0, q0, ctxt, n)

     alpha = 1.0_wp
     beta  = 0.0_wp

     do j = 1, n
       call pdgeadd ('notrans', n, 1, alpha, z, 1, &
            j, descz, beta, x, 1, 1, descx)
       call BLACS_BARRIER (ctxt, 'all')
       if (io_processor) then   
          write (nus) j
          write (nus) x
!          write(fo,*) x
!          vector_index=itoc(j)
!          len=lenth(vector_index)
!          Call IOsys('write real eigenvectors to eigen_file'//     &
!                     ' without rewinding',n,x,0,' ') 
       end if
     end do 
     if (io_processor) then        
!         Call IOsys('rewind eigenvectors on eigen_file read-and-write',0,0,0,' ')
!         Call IOsys('rewind all on eigen_file read-and-write',0,0,0,' ')
!         Call IOsys('close eigen_file',0,0,0,' ')
     end if
    deallocate (x)
    if (io_processor) endfile (unit=nus)
  end subroutine write_evec
