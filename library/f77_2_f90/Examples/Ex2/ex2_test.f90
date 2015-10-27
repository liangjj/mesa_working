  SUBROUTINE ex2_test(ig,n,inp,iout)
  INTEGER, DIMENSION(:), POINTER ::ig
  INTEGER n, iout, nloc
  read(inp,*) n
  write(iout,*) 'Reading in ',n,' as dimension'
  write(iout,*) 'Allocating ',n,' integers and filling array'
  allocate(ig(2*n))
  DO i=1,20
   ig(i)=i
  END DO 
  write(iout,*) 'Returning to calling routine and printing and deallocating'
  RETURN
  END SUBROUTINE ex2_test
















