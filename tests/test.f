      program test
      implicit none
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
      DO i=1,100
         a(i) = i
      END DO
      write(6,*) (a(i),i=1,100)
      stop
      end
