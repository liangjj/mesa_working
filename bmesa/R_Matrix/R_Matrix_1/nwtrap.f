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
      subroutine nwtrap(x1,x2,rtsafe,f,convg,niter,nactul,eigval,
     1                  gamma,n)
      implicit integer (a-z)
      real*8 x1, x2, convg, fl, df, fh, xl, xh, rtsafe, dxold
      real*8 dx, f, temp, rmat, drmat, eigval, gamma, kappa
      real*8 f1, f2
      dimension eigval(n), gamma(n)
      common/io/inp, iout
      rmat=0.d0
      call conrmat(rmat,rmat,eigval,gamma,gamma,x1,1,1,1,n,
     1            .false.,.false.)
      kappa=sqrt(-2.d0*x1)
      call fobj(fl,fl,rmat,rmat,kappa,1,.false.)
      rmat=0.d0
      call conrmat(rmat,rmat,eigval,gamma,gamma,x2,1,1,1,n,
     1             .false.,.false.)
      kappa=sqrt(-2.d0*x2)
      call fobj(fh,fh,rmat,rmat,kappa,1,.false.)
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
      f1=0.d0
      f2=0.d0
      call conrmat(f1,f2,eigval,gamma,gamma,rtsafe,1,1,1,n,
     1             .false.,.true.)
      kappa=sqrt(-2.d0*rtsafe)
      call fobj(f,df,f1,f2,kappa,1,.true.)
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
         f1=0.d0
         f2=0.d0
         call conrmat(f1,f2,eigval,gamma,gamma,rtsafe,1,1,1,n,
     1                .false.,.true.)
         kappa=sqrt(-2.d0*rtsafe)
         call fobj(f,df,f1,f2,kappa,1,.true.)
         if (f.lt.0.d0) then
             xl=rtsafe
         else
             xh=rtsafe
         endif
c         write(iout,3) iter, rtsafe, f
   10 continue
      write (iout,2) niter                                              
      return
    1 format (/,5x,'root not bracketed. will quit.',/,5x,
     1             'xl = ',e15.8,1x,'fl = ',e15.8/,5x,
     2             'xr = ',e15.8,1x,'fr = ',e15.8)
    2 format (/,5x,'no convergence after ',i4,' iterations') 
 3    format(/,1x,'iteration = ',i4,1x,'energy = ',e15.8,1x,
     1            'function = ',e15.8)
      end
