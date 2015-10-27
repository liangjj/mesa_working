!======================================================================
  module dc_matrix
!======================================================================

  use blacs, only: dlen_
  use precisn, only: wp

  implicit none

  real(wp), allocatable, save :: hmat(:) ! distributed Hamiltonian
  real(wp), allocatable, save :: z(:)    ! distributed Eigenvectors
  real(wp), allocatable, save :: b(:)    ! distributed Overlaps

  integer, save               :: desch(dlen_) ! Hamiltonian descriptor
  integer, save               :: descz(dlen_) ! Eigenvector descriptor
  integer, save               :: descb(dlen_) ! Overlap descriptor

  integer, save :: klsp

  real(wp), allocatable, save :: eval(:) ! eigenvalues

  end module dc_matrix
