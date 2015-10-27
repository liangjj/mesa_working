  SUBROUTINE rdmat
  USE lanczos_global
  IMPLICIT NONE 
  LOGICAL                                :: posinp
  INTEGER                                :: i,j
  if( posinp('$matrix',cpass)) then
      do i=1,n
         read(inp,*) (matrix(i,j), j=1,n)
      end do
  else
      stop
  end if
  title='input matrix'
  call prntcm(title,matrix,n,n,n,n,iout)
END SUBROUTINE rdmat
