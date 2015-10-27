*deck @(#)tolggr.f	1.1 9/8/91
      subroutine tolggr(funin,x,driver,funout,dfnout,l,scale,nbig,
     1                  subint,nsmall)
      implicit integer (a-z)
      real *8 driver, dfnout, scale, x
      complex *16 funin, funout
      common /io/ inp,iout
      dimension funin(nbig), funout(nsmall), driver(nbig)
      dimension dfnout(nsmall), x(nbig)
c----------------------------------------------------------------------c
c           put function and the effect of ( h  -e ) on function       c
c                                             0                        c
c                             on the small grid                        c
c----------------------------------------------------------------------c
      ipow=1
      if (l.eq.0) then
          ipow=0
      endif
      count=0
      do 10 i=1,nbig,subint
         count=count+1
         funout(count)=funin(i)
         dfnout(count)=driver(i)/scale
         dfnout(count)=dfnout(count)/(x(i)**ipow)
   10 continue
      return
      end

