*deck m6020
c***begin prologue     m6020
c***date written       910803   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m6020, link 6020, spline
c***author             rescigno, t. n.(llnl)
c***source             m6020
c***purpose            spline fit of free wave kohn functions
c***                   for zero potential. these are bessel type
c***                   solutions.
c***description        universal ( as function of rho=k*r) regular
c***                   and regularized irregular functions are computed,
c***                   and spline fit. spline coefficients put on file.
c***                   this version uses a much better algorithim for
c***                   computing the irregular function by solving an
c***                   inhomogeneous differential equation based on the
c***                   first born iterate for an exponential potential.
c***                   this produces a much smoother cutoff function which
c***                   still behaves properly at large rho.
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m6020
      program fitbes
      implicit integer (a-z)
      parameter (acc=30)
      character *8 cpass, bessfl, chrkey, filkne
      character *4096 ops
      character *800 card
      logical logkey, prnbes, prnspn, where
      real *8 z, fpkey, rmin, rmax, rdel, rd26, smldel, dl
      real *8 rl, alpha, jnint, ynint, jnintm, ynintm, scale
      real *8 energy, r0, r1, fuzz
      dimension z(1)
      common a(1)
      common /io/ inp, iout
      common /memory/ ioff
      equivalence(a,z)
      data energy, fuzz/ .5d+00,1.d-10 /
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      call posinp('$BESS',cpass)
      call cardin(card)
      lmax=intkey(card,'maximum-l-for-bessel',10,' ')
      nper=intkey(card,'points-per-interval-for-spline-coefficients',
     1                  4,' ')      
      ltop=intkey(card,'top-recursion-l',100,' ')       
      bessfl=chrkey(card,'bessel-file-name','BESSEL',' ')
      alpha=fpkey(card,'exponential-parameter',1.d+00,' ')
      subint=intkey(card,'integration-subinterval',3,' ')
      where=logkey(ops,'m6020=from-input',.false.,' ')
      if (.not.where) then
          call iosys ('read character "kohn data file" from rwf',0,
     1                 0,0,filkne)
          call iosys ('open kohndt as old',0,0,0,filkne)
          call iosys ('read real rmin from kohndt',1,rmin,0,' ')
          call iosys ('read real rmax from kohndt',1,rmax,0,' ')
      else
          rmin=fpkey(card,'rmin',0.d+00,' ')          
          rmax=fpkey(card,'rmax',10.d+00,' ')          
      endif
      prnbes=logkey(ops,'print=m6020=bessel',.false.,' ')
      prnspn=logkey(ops,'print=m6020=spline',.false.,' ')
      call iosys ('open bessel as new on ssd',262144,0,0,bessfl)
      call iosys ('write character "bessel function filename" to rwf',
     1             0,0,0,bessfl)
      call iosys ('write character "greens function" to rwf',0,0,0,
     1             'yes')
c---------------------------------------------------------------------c
c   note that rmin is the first spline knot point. if the real        c
c    grid goes below rmin, the cubic spline is being extrapolated     c
c                         (should be o.k.)                            c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                      generate spline data                           c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c             llmax = max l for splined bessel fcns                   c
c             rmax(rmin) = max(min) argument for                      c
c                          splined bessel fcns                        c
c             nper is the number of points to use per unit interval   c
c                         in computing the spline coefficients        c
c---------------------------------------------------------------------c  
      write (iout,100)
c---------------------------------------------------------------------c
c            estimate starting l                                      c
c---------------------------------------------------------------------c
      rl=lmax*acc
      strtl=lmax+sqrt(rl)
      strtl=max(strtl,lmax+1)
      ltop=min(ltop,strtl)
c-----------------------------------------------------------------------c
c       compute total number of points,nr,in spline                     c
c       compute the number of points in the integration.                c
c       this number is chosen such that the integration step size       c
c       will produce functional values at the spline grid points.       c       
c-----------------------------------------------------------------------c
      nr=nper*(rmax-rmin)
      nint=nr*subint
      rdel=(rmax-rmin)/nr
      smldel=(rmax-rmin)/nint
      nr=nr+1
      nint=nint+1
      rd26=rdel*rdel/6.d+00
      rmin=rmin+fuzz
      write(iout,10) rmin,rmax,lmax,ltop
      write (iout,20) nr, nint, rdel, smldel
      dimr=nr*(lmax+1)
      dimc=2*dimr
      words = 9*nint + 4*(ltop+1)*2 + 2 + 2*nr + 2*dimc +2*dimr
      write (iout,30) words
      call iosys ('read integer maxsiz from rwf',1,canget,0,' ')
      if (words.gt.canget) then
          call lnkerr ('cannot get required memory:will quit')
      endif
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(words,z,ngot,'m6020',0)
      x=ioff
      xinv=x+nint
      lastx=xinv-1
      expfn=xinv+nint
      lastxi=expfn-1
      v=expfn+nint
      j=v+nint
      jp=j+2*(ltop+1)
      y=jp+2*(ltop+1)
      yp=y+2*(ltop+1)
      norm=yp+2*(ltop+1)
      psi0=norm+2
      psi1=psi0+nint
      driver=psi1+nint
      fun=driver+nint
      fnl=fun+2*nint
      ddfnl=fnl+dimc
      coefc=ddfnl+dimr
      coefr=coefc+dimc
      scr=coefr+dimr
c----------------------------------------------------------------------c
c           make integration grid and exponential function             c
c----------------------------------------------------------------------c
      call mkxexp(z(x),z(xinv),z(expfn),alpha,rmin,smldel,nint)
c
      locfn=fnl
      locdfn=ddfnl
      do 200 l=0,lmax
c----------------------------------------------------------------------c
c              calculate effective potential for this l                c
c----------------------------------------------------------------------c
         call potntl(z(v),rmin,l,smldel,nint)
c----------------------------------------------------------------------c
c        integrate the homogeneous equation for the bessel             c
c                        function                                      c
c----------------------------------------------------------------------c
         call rzero(z(driver),nint)
         r0=rmin**(l+1)
         r1=(rmin+smldel)**(l+1)
         call numerv(z(psi0),r0,r1,z(driver),energy,z(v),smldel,
     1               nint)
c----------------------------------------------------------------------c
c              calculate regular and irregular bessel                  c
c              functions to normalize wavefunction                     c
c----------------------------------------------------------------------c
         dl=l
         if ( dl.gt.rmax) then
c----------------------------------------------------------------------c
c             backward recursion                                       c
c----------------------------------------------------------------------c
              call rcbesb(z(lastx-1),z(lastxi-1),z(j),z(jp),z(y),z(yp),
     1                    z(norm),2,2,l,ltop,.false.)
         else
c----------------------------------------------------------------------c
c              forward recursion                                       c
c----------------------------------------------------------------------c
              call rcbesf(z(lastx-1),z(lastxi-1),z(j),z(jp),z(y),
     1                    z(yp),2,2,l,ltop,.false.)
         endif
         jnint=z(j+2*l+1)
         ynint=z(y+2*l+1) 
         jnintm=z(j+2*l)
         ynintm=z(y+2*l)    
c----------------------------------------------------------------------c
c             calculate new driver                                     c
c----------------------------------------------------------------------c
         call drver(z(driver),z(expfn),z(psi0),nint)
c----------------------------------------------------------------------c
c               integrate inhomogeneous differential equation          c
c----------------------------------------------------------------------c
         call numerv(z(psi1),r0,r1,z(driver),energy,z(v),smldel,
     1               nint)
c----------------------------------------------------------------------c
c              compute the constants necessary to make                 c
c              the solution go to a pure neumann function.             c
c              then construct the final complex function.              c
c----------------------------------------------------------------------c
         call fnlfun(z(fun),scale,z(psi0),z(psi1),z(x),jnint,ynint,
     1               jnintm,ynintm,nint)
c----------------------------------------------------------------------c
c         put function and second derivative to be interpolated        c
c                        on to small grid                              c
c----------------------------------------------------------------------c
         call tolggr(z(fun),z(x),z(driver),z(locfn),z(locdfn),l,scale,
     1               nint,subint,nr)
c----------------------------------------------------------------------c
c                  up the counters                                     c
c----------------------------------------------------------------------c
         locfn=locfn+2*nr
         locdfn=locdfn+nr
  200 continue
c----------------------------------------------------------------------c
c               calculate x on small grid                              c
c----------------------------------------------------------------------c
      call mkx(z(x),rmin,rdel,nr)
c----------------------------------------------------------------------c
c                write out necessary stuff to bessel                   c
c----------------------------------------------------------------------c      
      call iosys ('write integer maximum-l to bessel',1,lmax,0,' ')
      call iosys ('write integer "total pts" to bessel',1,nr,0,' ')
      call iosys ('write real points to bessel',nr,z(x),0,' ')
      call iosys ('write real rmin to bessel',1,rmin,0,' ')
      call iosys ('write real rmax to bessel',1,rmax,0,' ')
      call iosys ('write real "r spacing" to bessel',1,rdel,0,' ')
c----------------------------------------------------------------------c
c              make the spline coefficients                            c
c----------------------------------------------------------------------c
      locfn=fnl
      locdfn=ddfnl
      loccmp=coefc
      locrel=coefr
      do 300 l=0,lmax      
         call splinc(z(x),z(locfn),rdel,nr,z(scr),z(loccmp))
         call splinr(z(x),z(locdfn),rdel,nr,z(scr),z(locrel))
         locfn=locfn+2*nr 
         locdfn=locdfn+nr
         loccmp=loccmp+2*nr
         locrel=locrel+nr
  300 continue
      call iosys('write real "complex kohn function" to bessel',
     1            dimc,z(fnl),0,' ')            
      call iosys('write real "second derivative of complex kohn '//
     1           'function" to bessel',dimr,z(ddfnl),0,' ')            
      call iosys('write real "function spline coefficients" to '//
     1           'bessel',dimc,z(coefc),0,' ')            
      call iosys('write real "second derivative spline coefficients" '//
     1            'to  bessel',dimr,z(coefr),0,' ')            
      if (prnspn) then
          call prspln(z(x),z(fnl),z(ddfnl),z(coefc),z(coefr),nr,lmax)
      endif
   10 format(/,5x,'rmin spline',1x,f10.5,1x,'rmax spline',1x,f10.5,1x,
     1            'l max',1x,i3,1x,'l top',1x,i3)
   20 format(/,5x,'no. spline points',1x,i5,1x,'no. integration points',
     1         1x,i5,/,5x,'spline spacing',1x,f10.5,
     2                 1x,'integration step',1x,f10.5)
   30 format (/,5x,'get',1x,i8,1x,'words of memory')
  100 format(//,20x,'***** m6020:splined ricatti-bessel functions *****'
     1)
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call iosys ('rewind all on bessel read-and-write',0,0,0,' ')
      call iosys ('close bessel',0,0,0,' ')
      call chainx(0)
      stop
      end
*deck tolggr
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
*deck fnlfun
c***begin prologue     fnlfun
c***date written       910805   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m6020, link 6020, spline
c***author             schneider b.(lanl)
c***source             m6020
c***purpose            construct complex kohn free wave
c***                   
c***description        a single, complex wave whose imaginary part
c***                   is a regular ricatti-bessel function for all
c***                   values of rho and whose real part is a ricatti-
c***                   neumann function asymptotically, is constructed from
c***                   a bessel function and a solution of the first born
c***                   iterate to an exponential potential. 
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       fnlfun
      subroutine fnlfun(final,scale,psil0,psil1,x,jnint,ynint,jnintm,
     1                  ynintm,n)
      implicit integer (a-z)
      real *8 wron, psil0, psil1, jnint, ynint, x
      real *8 jnintm, ynintm, scale
      complex *16 final, ai, f1, f2, c1, c2
      dimension psil0(n), psil1(n), final(n), x(n)
c----------------------------------------------------------------------c
c                 the solution needed is:                              c
c                    soln = (i*j  - q ) / x                            c
c                               l    l                                 c
c            where q behaves asymptoticall as a neumann function       c
c            but is regular at the origin                              c 
c----------------------------------------------------------------------c
      ai=cmplx(0.d+00,1.d+00)
      wron=psil0(n)*psil1(n-1)-psil0(n-1)*psil1(n)
      f1=( -ynint + ai*jnint )
      f2=( -ynintm + ai*jnintm )
      c1=( f1*psil1(n-1) - f2*psil1(n) )/wron
      c2=( f2*psil0(n) - f1*psil0(n-1) )/wron
      do 10 i=1,n
         final(i) = ( c1*psil0(i) + c2*psil1(i) ) / x(i)
   10 continue
      scale=1.d+00/c2
      return
      end
*deck splinc
c***begin prologue     splinc
c***date written       890404   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           spline, link 6020, orbital decomposition
c***author             rescigno, t. n.(llnl)
c***source             spline
c***purpose            spline fit
c***references       
c
c***routines called    none
c***end prologue       spline
      subroutine splinc(x,y,del,n,scr,y2)
      implicit integer (a-z)
      real *8 x, del
      complex*16 y, yp1, ypn, scr, y2, p, qn, un, sig
      dimension x(n), y(n), scr(n), y2(n)
      yp1= ( y(2) - y(1) )/del
      ypn= ( y(n) - y(n-1) )/del
      y2(1)=-.5d+00
      scr(1)=(3.d+00/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      do 10 i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.d+00
         y2(i)=(sig-1.)/p
         scr(i)=(6.d+00*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     1           /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*scr(i-1))/p
   10 continue
      qn=.5d+00
      un=(3.d+00/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      y2(n)=(un-qn*scr(n-1))/(qn*y2(n-1)+1.d+00)
      do 20 k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+scr(k)
   20 continue
      return
      end
*deck splinr
c***begin prologue     splinr
c***date written       890404   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           spline, link 6020, orbital decomposition
c***author             rescigno, t. n.(llnl)
c***source             spline
c***purpose            spline fit
c***references       
c
c***routines called    none
c***end prologue       spline
      subroutine splinr(x,y,del,n,scr,y2)
      implicit integer (a-z)
      real *8 x, del
      real *8 y, yp1, ypn, scr, y2, p, qn, un, sig
      dimension x(n), y(n), scr(n), y2(n)
      yp1= ( y(2) - y(1) )/del
      ypn= ( y(n) - y(n-1) )/del
      y2(1)=-.5d+00
      scr(1)=(3.d+00/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      do 10 i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.d+00
         y2(i)=(sig-1.)/p
         scr(i)=(6.d+00*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     1           /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*scr(i-1))/p
   10 continue
      qn=.5d+00
      un=(3.d+00/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      y2(n)=(un-qn*scr(n-1))/(qn*y2(n-1)+1.d+00)
      do 20 k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+scr(k)
   20 continue
      return
      end
      subroutine mkxexp(x,xinv,expfn,alpha,rmin,rdel,n)
      implicit integer (a-z)
      real *8 x, xinv, expfn, rmin, rdel, alpha, one
      dimension x(n), expfn(n), xinv(n)
      data one/1.d+00/
      x(1)=rmin
      xinv(1)=one/x(1)
      expfn(1)=exp(-alpha*rmin)
      do 10 i=2,n
         x(i)=x(i-1) + rdel
         xinv(i)=one/x(i)
         expfn(i)=exp(-alpha*x(i))
   10 continue
      return
      end
      subroutine mkx(x,rmin,rdel,n)
      implicit integer (a-z)
      real *8 x, rmin, rdel
      dimension x(n)
      x(1)=rmin
      do 10 i=2,n
         x(i)=x(i-1) + rdel
   10 continue
      return
      end
      subroutine prspln(x,hs,hsder,cj,cy,nr,lmax)
      implicit integer (a-z)
      common/ io/ inp, iout
      complex *16 hs, cj
      real *8 x, hsder, cy
      dimension x(nr), hs(nr,0:lmax), hsder(nr,0:lmax)
      dimension cj(nr,0:lmax), cy(nr,0:lmax)
      do 10 l=0,lmax
         write (iout,100) l
         write (iout,200) (hs(i,l), i=1,nr)
         write (iout,110) l
         write (iout,200) (hsder(i,l), i=1,nr)
         write (iout,120) l
         write (iout,200) (cj(i,l), i=1,nr)
         write (iout,130) l
         write (iout,200) (cy(i,l), i=1,nr)
   10 continue
  100 format(/,5x,'hs for l',1x,i3)
  110 format(/,5x,'hsder for l',1x,i3)
  120 format(/,5x,'cj for l',1x,i3)
  130 format(/,5x,'cy for l',1x,i3)
  200 format(/,5x,4d15.8)
      return
      end
*deck drver
      subroutine drver(driver,expfn,psi,n)
      implicit integer (a-z)
      real *8 driver, expfn, psi
      dimension driver(n), expfn(n), psi(n)
      do 10 i=1,n
         driver(i)=expfn(i)*psi(i)
   10 continue
      return
      end        
*deck numerv
c***begin prologue     numerv
c***date written       910803
c***revision date               (yymmdd)
c***keywords           numerov, integration
c***author             schneider, barry(lanl)
c***source             @(#)util
c***purpose            to integrate one dimensional second order
c***                   differential equation for a regular solution
c***                   using three point numerov scheme.
c***
c***references         numerov method is well known and can be found
c***                   in many texts on numerical anaysis.
c
c***routines called  
c***end prologue
      subroutine numerv(psi,r0,r1,driver,energy,pot,stp,n)
      implicit integer (a-z)
      real *8  psi, driver, energy, stp, pot
      real *8  f, g, r0,r1
      dimension psi(n), pot(n), driver(n)
      common /io/ inp, iout
      f(i) = 2.d+00*( energy - pot(i) )
      g(i) = -2.d+00*driver(i)
c----------------------------------------------------------------------c
c               start solution at origin as regular                    c
c               with second point chosen as 1.d+00                     c
c----------------------------------------------------------------------c
      psi(1)=r0
      psi(2)=r1
c----------------------------------------------------------------------c
c               continue using numerov algorithim                      c
c----------------------------------------------------------------------c
      do 10 i=3,n
         psi(i)= ( 24.d+00 -10.d+00*stp*stp*f(i-1) )*psi(i-1)
         psi(i) = psi(i) - ( 12.d+00 + stp*stp*f(i-2) )*psi(i-2)
         psi(i) = psi(i) + stp*stp*( g(i) + 10.d+00*g(i-1) + g(i-2) )
         psi(i) = psi(i)/ ( 12.d+00 + stp*stp*f(i) )
   10 continue
      return
      end
*deck potntl
c***begin prologue     potntl
c***date written       910128   (yymmdd)
c***revision date               (yymmdd)
c***keywords           potential
c***author             schneider, barry(lanl)
c***source             @(#)m6020
c***purpose            calculate potential for one dimensional
c***                   schroedinger equation.
c***
c***
c***references
c
c***routines called    util
c***end prologue
      subroutine potntl(v,r0,l,stp,n)
      implicit integer (a-z)
      dimension v(n)
      real *8 v, stp, rval, r0
      common /io/ inp, iout
      if (r0.eq.0.0d+00) then
          r0=1.d-06
      endif
      rval=r0
      do 10 i=1,n
         v(i)=.5d+00*l*(l+1)/(rval*rval)
         rval=rval+stp
   10 continue
      return
      end
