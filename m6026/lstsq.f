*deck lstsq
c***begin prologue     lstsq
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
c                      pt are the points
c                      rhs is the output set of coefficients
c                      the rest of the variables are used only internally
c***references
c***routines called    sgemm(clams)
c***end prologue       lstsq
      subroutine lstsq(f,pt,rhs,xn,coef,ipvt,npnts,nlpar,rms)
      implicit integer (a-z)
      real *8 f, coef, rhs, xn, pt, rms, sdot, value
      character *80 title
      dimension f(npnts), pt(npnts), rhs(nlpar)
      dimension  xn(nlpar,npnts), coef(nlpar,nlpar), ipvt(nlpar)
      common /io/ inp, iout
c----------------------------------------------------------------------c
c                     print out input information                      c
c----------------------------------------------------------------------c
      write (iout,100) nlpar
  100 format(/,5x,'least squares polynomial fitting subroutine',//,5x,
     1            'order of fit',1x,i3)
      call rzero(xn,nlpar*npnts)
c----------------------------------------------------------------------c
c                    set up coefficient matrix                         c
c----------------------------------------------------------------------c
      do 10 j=1,npnts
         xn(1,j)=1.d0
   10 continue
      if (nlpar.gt.1) then
          do 20 i=2,nlpar
             do 30 j=1,npnts
                xn(i,j)=xn(i-1,j)*pt(j)
   30        continue
   20     continue      
      endif
      call rzero(coef,nlpar*nlpar)
      call ebct(coef,xn,xn,nlpar,npnts,nlpar)
c----------------------------------------------------------------------c
c               factor coefficient matrix                              c
c----------------------------------------------------------------------c
      call sgefa(coef,nlpar,nlpar,ipvt,info)
      if (info.ne.0) then
          write(iout,1)
      endif
c----------------------------------------------------------------------c
c                get right hand side                                   c
c----------------------------------------------------------------------c
      call ebc(rhs,xn,f,nlpar,npnts,1)
c----------------------------------------------------------------------c
c               solve equations for coefficients                       c
c----------------------------------------------------------------------c
      call sgesl(coef,nlpar,nlpar,ipvt,rhs,0)
c----------------------------------------------------------------------c
c            calculate the rms error                                   c
c----------------------------------------------------------------------c
      rms=0.d0
      do 40 i=1,npnts
         value=sdot(nlpar,rhs,1,xn(1,i),1)
         rms=rms+(f(i)-value)*(f(i)-value)
   40 continue
      write(iout,*) 'rms error',rms
      title='least squares coefficients'
      call prntrm(title,rhs,nlpar,1,nlpar,1,iout)
      return
    1 format(/,5x,'singular matrix in least squares routine')
      end

