*deck @(#)fmch.f	1.1 9/7/91
c***begin prologue     fmch
c***date written                (yymmdd)
c***revision date      890417   (yymmdd)
c***keywords           fmch, link 6003, error function
c***author             unknown
c***source             m6002
c***purpose            special functions for gaussian integrals
c*** 
c
c***references       
c
c***routines called    
c***end prologue       fmch
      function fmch(m,x,y)
      implicit real *8 (a-h,o-z)
      common /io/ inp, iout
      if (x-10.d0) 10,10,20
   10 a=m
      a=a+0.5d0
      term=1.0d0/a
      ptlsum=term
      do 11 i=2,50
         a=a+1.0d0
         term=term*x/a
         ptlsum=ptlsum+term
         if (term/ptlsum-0.00000001d0) 12, 11, 11
   11 continue
      write (iout,999) m,x
      call lnkerr('error in fmch')
   12 fmch=0.5d0*ptlsum*y
      go to 150
   20 a=m
      b=a+0.5d0
      a=a-0.5d0
      xd=1.d0/x
      approx=0.88622692d0*( sqrt(xd)*xd**m)
      if (m) 21, 23, 21
   21 do 22 i=1,m
         b=b-1.0d0
         approx=approx*b
   22 continue
   23 fimult=0.5d0*y*xd
      ptlsum=0.d0
      if (fimult) 421,25,421
  421 continue
      fiprop=fimult/approx
      term=1.0d0
      ptlsum=term
      notrms=x
      notrms=notrms+m
      do 24 i=2,notrms
         term=term*a*xd
         ptlsum=ptlsum+term
         if ( (abs(term*fiprop/ptlsum)-0.00000001d0).gt.0.d+00) then
              a=a-1.0d0
         else
              go to 25
         endif
   24 continue
      write (iout,999) m,x
      call lnkerr('error in fmch')      
   25 fmch=approx-fimult*ptlsum
  150 return
  999 format (/,5x,'no convergence for fmch', i6, e16.9)
      end
