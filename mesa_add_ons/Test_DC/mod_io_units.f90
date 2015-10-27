!----------------------------------------------------------------------
      Module io_units
!----------------------------------------------------------------------
! ... define unit numbers for standard input and output
!----------------------------------------------------------------------

      implicit none

      integer, parameter      :: fi = 8    ! standard input
      character(40)           :: AF_inp = 'pgd.inp'
      integer, parameter      :: fo = 9    ! standard output
      character(40)           :: AF_out = 'pgd.out'
      integer, parameter      :: nuh = 80   ! matrix input
      character(40)           :: AF_int = 'dc_mat.nnn'
      integer, parameter      :: nus = 90   ! solution output
      character(40)           :: AF_sol = 'dc_sol.nnn'
      character(LEN=240)      :: card
      character(LEN=80)       :: cpass
      End module io_units
