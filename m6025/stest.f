      subroutine stest ( x, a, n, s, iflag, avg, sdev )
      implicit real *8 (a-h,o-z)
c
c      tests to see if x is within s standard deviations
c       of the average in the array a
c
      dimension a(n)
c
      sum = 0d0
      sum2 = 0d0
c
      do 1 i = 1, n
         sum = sum + a(i)
         sum2 = sum2 + a(i)**2
 1    continue
c
      avg = sum/n
      dev = sum2/n
      sdev = sqrt(float(n)*(dev - avg*avg)/(n-1))
c
      iflag = 1
      if ( abs(x-avg) .lt. s*sdev ) iflag = 0
c
c      print 2000, avg, sdev, iflag
c 2000 format ( " result = avg,sdev,iflag ", 2g12.4, i6 )
      return
      end
