*deck nwtrap
c***begin prologue     nwtrap
c***date written       930623   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           newton-raphson root finder
c***author             schneider, barry (nsf)
c***source             math
c***purpose            find zero of function using combination of bisection
c***                   and newton-raphson method
c***                   
c***                   
c
c***references         
c
c***routines called    
c***end prologue       nwtrap
      subroutine nwtrap(x1,x2,rtsafe,convg,niter,nactul,f1,f2,f3,f4,
     1                  l,ltop,type)
      implicit integer (a-z)
      real*8 x1, x2, convg, fl, df, fh, xl, xh, rtsafe, dxold
      real*8 dx, f, temp, f1, f2, f3, f4
      character*(*) type
      dimension f1(0:ltop), f2(0:ltop), f3(0:ltop), f4(0:ltop)
      common/io/inp, iout
      call rc1bes(x1,f1,f2,f3,f4,l,ltop,'derivatives',.false.)
      fl=f1(l)
      if(type.eq.'jp') then
         fl=f2(l)
      endif
      call rc1bes(x2,f1,f2,f3,f4,l,ltop,'derivatives',.false.)
      fh=f1(l)
       if(type.eq.'jp') then
         fh=f2(l)
      endif
      if (fl*fh.ge.0.d0) then
          write (iout,1) x1, fl, x2, fh
          call lnkerr('quit')
      endif
      if (fl.lt.0.d0) then
           xl=x1
           xh=x2
      else
           xl=x2
           xh=x1
      endif
      rtsafe=.5d0*(x1+x2)
      dxold=abs(x2-x1)
      dx=dxold
      call rc1bes(rtsafe,f1,f2,f3,f4,l,ltop,'derivatives',.false.)
      f=f1(l)
      df=f2(l)
      if (type.eq.'jp') then
          f=f2(l)
          call secder(f1(l),df,rtsafe,l,1)
      endif
      do 10 iter=1,niter
         if (((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0.d0.or.
     1         abs(2.d0*f).gt.abs(dxold*df)) then
               dxold=dx
               dx=.5d0*(xh-xl)
               rtsafe=xl+dx
               if ( xl.eq.rtsafe ) then
                    nactul=iter
                    return
               endif     
         else
               dxold=dx
               dx=f/df
               temp=rtsafe
               rtsafe=rtsafe-dx
               if ( temp.eq.rtsafe ) then
                    nactul=iter
                    return
               endif
         endif     
         if (abs(dx).lt.convg) then
             nactul=iter
             return
         endif
         call rc1bes(rtsafe,f1,f2,f3,f4,l,ltop,'derivatives',.false.)
         f=f1(l)
         df=f2(l)
         if (type.eq.'jp') then
             f=f2(l)
             call secder(f1(l),df,rtsafe,l,1)
         endif 
         if (f.lt.0.d0) then
             xl=rtsafe
         else
             xh=rtsafe
         endif
   10 continue
      write (iout,2) niter                                              
      return
    1 format (/,5x,'root not bracketed. will quit.',/,5x,
     1             'xl = ',e15.8,1x,'fl = ',e15.8/,5x,
     2             'xr = ',e15.8,1x,'fr = ',e15.8)
    2 format (/,5x,'no convergence after ',i4,' iterations') 
      end
