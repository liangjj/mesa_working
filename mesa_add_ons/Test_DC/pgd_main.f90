!======================================================================
  program pgd
!======================================================================
! Solution of Generalized Eigenproblem for DC matrices
!----------------------------------------------------------------------

  use precisn, only: wp
  use io_units
  use blacs, only: nprocs, iam, ictxt, io_processor, p, q, ctxt, myrow,&
       mycol, p0, q0, set_cur_ctxt, p_error, nblock
  use dc_matrix, only: klsp

  implicit none

  real(wp)                    :: t0, t1
  integer                     :: imycol, imyrow, ip, iq
  integer                     :: c0, c1, cr
  integer                     :: parms(4)
  integer                     :: intkey
  integer                     :: len
  integer                     :: lenth
  logical                     :: dollar
  character (len=80)          :: chrkey 
  character (len=3)           :: itoc 
  character (len=3)           :: suffix 

  call cpu_time (t0)
  call system_clock (count=c0)
!----------------------------------------------------------------------
! define liniar blacs grid to include all processes to broadcast:

  call BLACS_PINFO (iam, nprocs) ! find process #, total # processors
  call BLACS_GET (-1, 0, ictxt)  ! find default context, ictxt

  call BLACS_GRIDINIT (ictxt, 'Row-major', 1, nprocs)
  call BLACS_GRIDINFO (ictxt, ip, iq, imyrow, imycol)
  io_processor = (imycol == 0)
  call set_cur_ctxt (1)
  if (io_processor) then 
!     OPEN(fi,file=AF_inp,status='old')
     open (fi,file=AF_inp)
!     IF ( dollar('$general_eigen_solve',card,cpass,fi) ) then
!          p = intkey(card,'number_of_processors',2,' ')
!          q = intkey(card,'number_of_cores_per_processor',2,' ')
!          nblock = intkey(card,'number_of_matrix_blocks',64,' ')
!          klsp = intkey(card,'file_suffix',1,' ')
!          AF_out = chrkey(card,'output_file_name','pgd.out',' ')
!          len = lenth(AF_out)
!          OPEN(fo,file=AF_out(1:len),status='unknown')
          OPEN(fo,file=AF_out,status='unknown')
          write (fo,*)
          write (fo,'(15x,a)') '============='
          write (fo,'(15x,a)') ' Program PGD '
          write (fo,'(15x,a)') '============='
          write (fo,*)
          write (fo,'(a,i6)') 'Number of processors = ', nprocs
!     END IF
     Call Read_ipar(fi,'p',p)
     Call Read_ipar(fi,'q',q)
     Call Read_ipar(fi,'nblock',nblock)
     write (fo,*)
     write (fo,'(a,3i6)') 'Blacs process dimensions, p,q,nblock = ', &
                           p,q,nblock
     write (fo,*)
     Call Read_ipar(fi,'klsp',klsp)
     write (fo,'(a,i3)') 'DC Parallel General Eigensolution for klsp =',klsp
     parms = (/nblock, p, q, klsp/)
     call igebs2d (ictxt, 'all', ' ', 4, 1, parms, 4)

  else
     call igebr2d (ictxt, 'all', ' ', 5, 1, parms, 5, 0, 0)
     nblock = parms(1)
     p = parms(2)
     q = parms(3)
     klsp = parms(4)
  end if
  call BLACS_BARRIER (ictxt, 'All')
  call BLACS_GRIDEXIT (ictxt)     ! kill initial mesh

  if(p*q > nprocs) Stop 'p*q > nprocs'
!----------------------------------------------------------------------
! create p-q Blacs grid:

  call BLACS_GET (-1,0,ctxt)  ! find default context, ctxt
  call BLACS_GRIDINIT (ctxt, 'Row-major', p, q)
  call BLACS_GRIDINFO (ctxt, ip, iq, myrow, mycol)

  if (myrow >= 0 .and. myrow < p .and. mycol >= 0 .and. mycol < q) then

     io_processor = (myrow == p0 .and. mycol == q0)
     call set_cur_ctxt (2)

     call pgdiag             ! perform eigensolution

     if (io_processor) then
        write (fo,'(/,a,/)') 'end of PGD'
        call cpu_time (t1)
        write (fo,'(a,f16.4,a)') 'CPU time     = ', t1 - t0, ' secs'
        call system_clock (count=c1, count_rate=cr)
        write (fo,'(a,f16.4,a)') 'Elapsed time = ', REAL(c1-c0,wp) / &
             REAL(cr,wp), ' secs'
     end if

  end if

  call BLACS_BARRIER (ctxt, 'All')
  call BLACS_GRIDEXIT (ctxt)    
  call BLACS_EXIT(0)

  stop
  end program pgd
