module precisn
! define working precision kind value
! Time-stamp: "2001-03-20 09:46:41 cjn"
  implicit none
  public
  integer, parameter    :: wp = selected_real_kind(12)
end module precisn
