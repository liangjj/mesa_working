module io_units
! define unit numbers for standard input and output
! Time-stamp: "2008-10-23 15:42:37 cjn"
  implicit none
  private
  public fi, fo, nfnml, nuh
  integer, parameter      :: fi = 5    ! standard input
  integer, parameter      :: fo = 6    ! standard output
  integer, parameter      :: nfnml = 7 ! namelist input
  integer, parameter      :: nuh = 8   ! eigenvecotr output
end module io_units
