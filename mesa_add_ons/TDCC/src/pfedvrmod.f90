module pfedvrmod
  USE nrtype
  implicit NONE

  !----------------------ALLOCATABLE type definition---------------------------
  ! The dimensions here are (1:num_func,1:num_func) for the matrices
  ! and (ndim,1:num_func) for the vectors

  TYPE dvr_mat
     REAL(DP), dimension(:,:), pointer :: ke_mat, eigvec_mat, df, ddf
     REAL(DP), dimension(:), pointer   :: pt, wt, fac1, fac2, eigval_mat
  END TYPE dvr_mat


  !--------The length of mat_reg (a dvr_mat type) is (ndim)x(the number of regions in each dimension) ------------

  TYPE (dvr_mat), dimension(:,:), allocatable :: mat_reg 


  !------Number of dvr functions in each region & dimension------------------

  INTEGER(I4B), dimension(:,:), allocatable :: num_fun


  !------- Lower and upper boundaries of each dimension & region------------------
  !------- in the format (ndim,region,2)--> 2=x_min,x_max------------------------------

  REAL(DP), dimension(:,:,:), allocatable :: bounds  


  !-------------Number of dvr regions for each dimension--------------------------------------

  INTEGER(I4B), dimension(:), allocatable :: num_reg

end module pfedvrmod
