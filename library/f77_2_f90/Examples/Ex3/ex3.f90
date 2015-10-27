  PROGRAM test
  USE test_global
  IMPLICIT NONE
  INTEGER                               :: i, j, k
  open (inp,file='inp',status='old')
  open (iout,file='out',status='unknown')  
  read(inp,*) m, ( nr(i),i=1,m ), (nc(i),i=1,m)
  allocate(pg(m))
  offset=0
  do i=1,m
   allocate(pg(i)%g(nr(i),nc(i)),pg(i)%h(nr(i)))
   call ex3_test(pg(i)%g,pg(i)%h,nc(i),nr(i))
   write(iout,*) 'The ',i,' Matrix'
   write(iout,*) 'nr = ',nr(i),'nc = ',nc(i)
   do j=1,nr(i)
     write(iout,*) (pg(i)%g(j,k),k=1,nc(i))
   end do
   write(iout,*) 'The ',i,' Vector'
   write(iout,*) (pg(i)%h(k),k=1,nr(i))
   offset=offset+3
  end do
  call mmat
  stop
  deallocate(pg)
  END PROGRAM test









