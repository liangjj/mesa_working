program test
      implicit none
      integer, dimension(:), allocatable    :: itest
      integer maxmem, stat, err, i
      integer a, realloc
      pointer(p,a(1))
      maxmem=500000000
      p=malloc(4*maxmem)
      IF (p.eq.0) THEN
          write(6,*) 'Malloc Failure'
          stop
      END IF
      maxmem=maxmem+10000
      p=realloc(p,maxmem)
      ALLOCATE(itest(maxmem),stat=err)
      DO i=10000,11000
         itest(i) = i
      END DO
      write(6,*) itest(10000:11000)
      IF (err /= 0 ) THEN
          write(6,*) 'There is an allocate error'
          stop
      ELSE
          write(6,*) 'All is Well'
      END IF
      DEALLOCATE(itest)
end program test
