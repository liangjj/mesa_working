!deck fd_input
!***begin prologue     fd_input
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           finite difference, input
!***author             schneider, b. i.(nsf)
!***source             fdlib
!***purpose            input subroutine for finite difference library
!***references
!***routines called    
!***end prologue       fd_input

  SUBROUTINE fd_input(nphy,nglobal,row,coord)
  USE fd_global
  USE fd_prnt
  IMPLICIT NONE
  INTEGER                          :: nphy, nglobal, row
  CHARACTER(LEN=*)                 :: coord
  CHARACTER(LEN=80)                :: chrkey
  INTEGER                          :: intkey, lenth, len
  INTEGER                          :: i
  LOGICAL                          :: dollar, logkey
  REAL*8                           :: fpkey
  WRITE(iout,1) coord
  LEN=lenth(coord)
!
!     read the input by stopping at the keyword in the input stream.
!
  keywrd='$h0('//coord(1:LEN)//')'
  IF ( dollar(keywrd,card,cpass,inp) ) THEN
!
!      set the print variables for printing.
!
      prn_fd_loc = 'print='//prn_fd
      prn_fd_loc(7)=chrkey(card,'print',prn_fd(7),' ')
      IF(prn_fd_loc(7) == 'all') THEN
         prn_fd_log = .true.
      ELSE
         CALL setlog(prn_fd_log,prn_fd_loc,card,6)
      END IF
      units=chrkey(card,'units','atomic_units',' ')
      if(units == 'atomic_units') then
         hbar = 1.d0
         mass = 1.d0
      endif
      nodiag=logkey(card,'do_not_diagonalize',.false.,' ')
      ndiff=intkey(card,'order_of_finite_difference_formula',3,' ')
      row = ndiff/2 + 1
      CALL fparr(card,'region_boundaries',edge,2,' ')
!
!          Read in number of points
!  
      nstep=intkey(card,'number_of_steps',3,' ')
      del=(edge(2)-edge(1))/nstep
      angmom=intkey(card,'angular_momentum',0,' ')
      charge=fpkey(card,'charge',1.d0,' ')
  
!          The actual number of unspecified points is nphy.  The value of
!          the solution at the first and last point are fixed to be zero.
  
      nphy=nstep-1
      nglobal=nstep+1
      WRITE(iout,2) ndiff, nglobal, nphy, del
      WRITE(iout,3) (edge(i),i=1,2)
      write(iout,4) angmom, charge
  ELSE
      write(iout,5) keywrd
      stop
  END IF
1    FORMAT(/,20X,'coordinate = ',a24)
2    FORMAT(/,5x,'order of finite difference formula     = ',i2,  &
            /,5x,'number of points(including end points) = ',i5,  &
            /,5x,'number of physical points              = ',i5,  &
            /,5X,'step size                              = ',e15.8)
3    FORMAT(/,5X,'region edges = ',/,2(e15.5,1x))
4    FORMAT(/,5X,'angular momentum = ',i3, &
            /,5x,'charge           = ',e15.8) 
5    FORMAT(/,20x,'no input found corresponding to keyword = ',a80)
END SUBROUTINE fd_input
