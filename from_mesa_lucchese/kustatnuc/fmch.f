      function fmch(m,x,y)
      implicit real*8 (a-h,o-z)
      if (x-10.e0) 10,10,20
   10 a=m
      a=a+0.5e0
      term=1.0e0/a
      ptlsum=term
      do 11 i=2,50
      a=a+1.0e0
      term=term*x/a
      ptlsum=ptlsum+term
      if (term/ptlsum-0.00000001e0) 12, 11, 11
   11 continue
      write(66,999) m,x
      call exit
   12 fmch=0.5e0*ptlsum*y
      go to 150
   20 a=m
      b=a+0.5e0
      a=a-0.5e0
      xd=1.e0/x
      approx=0.88622692e0*( sqrt(xd)*xd**m)
      if (m) 21, 23, 21
   21 do 22 i=1,m
      b=b-1.0e0
   22 approx=approx*b
   23 fimult=0.5e0*y*xd
      ptlsum=0.e0
      if (fimult) 421,25,421
  421 continue
      fiprop=fimult/approx
      term=1.0e0
      ptlsum=term
      notrms=x
      notrms=notrms+m
      do 24 i=2,notrms
      term=term*a*xd
      ptlsum=ptlsum+term
      if ( abs(term*fiprop/ptlsum)-0.00000001e0)  25, 25, 24
   24 a=a-1.0e0
      write(66,999) m,x
      call exit
   25 fmch=approx-fimult*ptlsum
  150 return
  999 format (24h no convergence for fmch, i6, e16.9)
      end
