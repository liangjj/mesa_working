  SUBROUTINE mmat
 USE test_global
 write(iout,*) (nr(i),i=1,m)
  write(iout,*)( nc(i),i=1,m)
do i=1,m
  do j=1,nr(i)
     write(iout,*) ( pg(i)%g(j,k), k=1,nc(i))
  end do
end do
END SUBROUTINE mmat
















