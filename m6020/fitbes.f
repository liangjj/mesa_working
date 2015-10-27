*deck @(#)fitbes.f	1.1 9/8/91
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
      character *8 cpass, chrkey
      character *128 filkne, filbes
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
      call posinp('$bess',cpass)
      call cardin(card)
      lmax=intkey(card,'maximum-l-for-bessel',10,' ')
      nper=intkey(card,'points-per-interval-for-spline-coefficients',
     1                  4,' ')      
      ltop=intkey(card,'top-recursion-l',100,' ')       
      alpha=fpkey(card,'exponential-parameter',1.d+00,' ')
      subint=intkey(card,'integration-subinterval',3,' ')
      where=logkey(ops,'m6020=from-input',.false.,' ')
      if (.not.where) then
          call iosys ('read character "kohn data filename" from rwf',0,
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
      call iosys ('read character "bessel function filename" from rwf',
     1             -1,0,0,filbes)
      call iosys ('open bessel as new on ssd',262144,0,0,filbes)
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
      if (prnbes) then
          call prbes(z(x),z(fnl),z(ddfnl),nr,lmax)
      endif
      if (prnspn) then
          call prspln(z(x),z(coefc),z(coefr),nr,lmax)
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
