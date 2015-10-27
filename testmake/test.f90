  program test
  INTEGER                                  :: spdim, i, j, nreg    
  TYPE dvr_mat
    REAL*8, DIMENSION(:,:),   POINTER      :: ke_mat
  END TYPE dvr_mat
  TYPE (dvr_mat), DIMENSION(:,:), ALLOCATABLE &
                                           :: mat_reg
  spdim=3
  nreg=5
  ALLOCATE(mat_reg(nreg,i))
  DO i=1,spdim
  END DO
  DEALLOCATE(mat_reg)
  call exit
  end
