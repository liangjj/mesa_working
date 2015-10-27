subroutine propagate
  use globaali, only : method
  use globaali, only : rank              ! ### MPI ### 
  implicit none
  real*8 :: tbeg, tend
  ! ===================================================
  ! propagate the initial state in real/imaginary time
  ! ==================================================

  if(method .EQ. 'SO2') then
     call prop_so2
  else if(method .EQ. 'SO4') then
     call prop_so4
  endif
  
end subroutine propagate


