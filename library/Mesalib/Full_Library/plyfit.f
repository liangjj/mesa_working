*deck plyfit
c***begin prologue     plyfit
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           least squares,
c***author             schneider, barry (lanl)
c***source             mylib
c***purpose
c                      least squares polynomial fit of input function
c***description
c                      input function f is fit in a least squares sense
c                      to a power series
c
c                      f is input function
c                      power are the powers used in fitting
c                      pt are the points
c                      rhs is the output set of coefficients
c                      the rest of the variables are used only internally
c***references
c***routines called    sgemm(clams)
c***end prologue       plyfit
      subroutine plyfit(f,coef,rhs,xn,pt,power,ipvt,npwr,npnts,prnt,
     1                  first)
      implicit integer (a-z)
      real *8 f, coef, rhs, xn, pt, rms, sum
      logical prnt
      character *80 title
      character *(*) first
      dimension f(npnts), power(npwr), xn(npwr,npnts), pt(npnts)
      dimension coef(npwr,npwr), rhs(npwr), ipvt(npwr)
      common /io/ inp, iout
      if (first.eq.'first') then
c----------------------------------------------------------------------c
c                     print out input information                      c
c----------------------------------------------------------------------c
          write (iout,100) npwr
          write (iout,200) (power(i),i=1,npwr)
          call rzero(xn,npwr*npnts)
c----------------------------------------------------------------------c
c                    set up coefficient matrix                         c
c----------------------------------------------------------------------c
          do 10 i=1,npwr
             do 20 j=1,npnts
                xn(i,j)=pt(j)**power(i)
   20        continue
   10     continue
          if (prnt) then
              title='xn matrix'
              call prntrm(title,xn,npwr,npnts,npwr,npnts,iout)
          endif          
          call rzero(coef,npwr*npwr)
          call ebct(coef,xn,xn,npwr,npnts,npwr)
          if (prnt) then
              title='coefficient matrix'
              call prntrm(title,coef,npwr,npwr,npwr,npwr,iout)
          endif    
c----------------------------------------------------------------------c
c               factor coefficient matrix                              c
c----------------------------------------------------------------------c
          call sgefa(coef,npwr,npwr,ipvt,info)
          if (info.ne.0) then
              call lnkerr('error in linear system solve')
          endif    
      else
c----------------------------------------------------------------------c
c                get right hand side                                   c
c----------------------------------------------------------------------c
          call ebc(rhs,xn,f,npwr,npnts,1)
c----------------------------------------------------------------------c
c               solve equations for coefficients                       c
c----------------------------------------------------------------------c
          call sgesl(coef,npwr,npwr,ipvt,rhs,0)
          if (prnt) then
              title='least squares coefficients'
              call prntrm(title,rhs,npwr,1,npwr,1,iout)
          endif
c----------------------------------------------------------------------c
c              calculate rms error                                     c
c----------------------------------------------------------------------c
          rms=0.d0
          do 30 i=1,npnts
             sum=0.d0
             do 40 j=1,npwr
                sum=sum+rhs(j)*xn(j,i)
   40        continue
             rms=rms+f(i)-sum
   30     continue
          rms=sqrt(rms*rms)
          write (iout,300) rms                                       
      endif
      return
  100 format(/,5x,'least squares polynomial fitting subroutine',//,5x,
     1            'order of fit',1x,i3)
  200 format (//,5x,'powers',(/,15x,10(i2,1x)))
  300 format(/,1x,'rms error = ',e15.8)     
      end
