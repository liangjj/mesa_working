*deck resfit
      program resfit      
      implicit real *8 (a-h,o-z)
c
      character *8 ityp(2), name(3), ichan(3)
      character *10 rsn
      integer chan,type,initial,final
      logical fit, debug
      common /doyle/ npad, nit, apad(4), npad1(7), apad1(8), npad2(2)
      common /kar/ apad2(150), bl(30), bu(30), v(30), apad3(60), 
     1             nuflag(30), nbd(30)
      common / kchg  / apad4(300,15), bsav(300), nersdm
      common /kvar/ npar, nprob, ners, npad3(2), apad5, solves,
     1               npad4(2), apad6(2)
      common / pdata / np, npbck, elist(256), plist(256), sreal(256),
     >                 simag(256), phase(256), pabs(256), pres(256),
     >                 pbck(256), edup(256), ilist(256), del(256)
c
      namelist / res / emin,emax, chan,type,initial,final,fit,debug,
     >                 nit,npbck,escal,eref,rsn,elist,del,name,
     >                 np, xnpi
c
c     chan=1,2,3 for alpha,beta,gamma channels
c     type=1,2 for inelastic or reactive scattering
c     initial=0,1,2,... initial state label
c     final=0,1,2,...   final   state label
c     final = -1  ==>>  sum over all final states
c     emin,emax = range for spline fit
c        if not specified, then input table limits are used
c     escal = scale factor for energy
c     eref  = reference energy for spline fit
c             enew = escal*(eold - eref)
c     ebar  = list of special energies for output list
c     nbar  = number of elemnts in ebar list
c
c
c
      data ityp / 'inelastc', 'reactive' /
      data ichan / 'alpha', 'beta' , 'gamma' /
c
      open(unit=5,status='old',file='inp.res')
      open(unit=6,status='unknown',file='out.res')
      open(unit=9,status='unknown',file='scr.res')
c
c            establish default values
      emin = -999d0
      emax = -999d0
      pi = acos(-1d0)
      chan = 1
      type = 2
      initial = 0
      final = -1
      fit = .true.
      debug = .false.
      de = 0d0
      ne = -1
      escal = 1d0
      eref = 0d0
      nbar = 0
      npar = 3
      npbck = 2
      xnpi = 2.d0
      np = 0
      name(1) = ' '
      name(2) = ' '
      name(3) = ' '
      nit = 10
      rsn = ' '
c            read data
      nread = 0
 1    continue
      read (5,res,end=100)
c
      do 3 i = 1, npar
         read(5,*) nuflag(i), v(i), nbd(i),bl(i),bu(i)
 3    continue
      ntpar = npar + npbck + 1
c

      idup = 0
      do 52 i = 1, np
         sreal(i) = 0d0
         simag(i) = 0d0
         plist(i) = 0d0
         phase(i) = del(i)/pi
 52   continue
c
 51   continue
      call phasit ( elist, phase, np, xnpi, 0.d0,-1.d0, pres, pabs )
      do 11 i = 1, np
         pabs(i) = pabs(i)*pi
 11   continue
c
      write(6,1300) name,ityp(type),ichan(chan),initial,final,rsn
 1300 format ( '1sorted list of s-matrix elements for the transition',/,
     >         10x, 3a2,1x,a8,1x,a5, ' from v=',i2,'    to v='i2,
     >         ',   rsn=',a10, /,
     >          14x, 'energy     smat**2     s(real)     s(imag)',
     >          '       phase/pi    abs phase' )
      write(6,1301) (elist(i),plist(i),sreal(i),simag(i),phase(i),
     >               pabs(i),i=1,np)
 1301 format ( 5x, 0pf15.7, 1p3e12.3, 0pf15.4, f13.4 )
      if ( idup .gt. 0 ) print 1302, idup, (edup(i),i=1,idup)
 1302 format ( ' a total of ', i3,' duplicate energies found -- ',/,
     >         (1x,10f13.6) )
c
      write(6,1010) npbck, nit,
     >            (nuflag(i), v(i), nbd(i), bl(i), bu(i), i=1, 3)
 1010 format ( '1   non-linear least squares resonance fitting program',
     >         /,'0resonant part of phase shift fitted to 3 nonlinear'
     >           ' parameters',
     >         /,'      pres(e)= p*(.5*pi + atan(2*(e-eres)/gres))',
     >         /,'          parameters are -- p,eres,gres',
     >         /,'0background part of phase shift fitted to e**',i2,
     >         /,'0maximum number of iterations = ', i3,/,
     >         /,'0',12x,'vary flag    initial guess    bound flag',
     >           '    lower bound    upper bound',
     >         /,'0   -p-      ', i6, g20.6, i10, g19.6, g15.6,
     >         /,'   -eres-    ', i6, g20.6, i10, g19.6, g15.6,
     >         /,'   -gres-    ', i6, g20.6, i10, g19.6, g15.6 )
c
      ners = np
      init=0
      call opac (init,6,9,debug)
      call fcn ( v, bsav, phi )
c
      write(6,1011) (v(i),i=1,ntpar)
 1011 format ( '1final fit of data follows --', /,
     >        '0resonant amplitude      p    = ', g15.5, /,
     >        '0resonance energy     eres    = ', g15.5, /,
     >        '0resonance width      gres    = ', g15.5, /,
     >        '0linear background parameters = ', 5g15.5,/,(30x,5g15.5))
      write(6,1012)
 1012 format ( //,'0     energy    actual phase    fit phase    ',
     >         'background     resonance         error' )
      do 15 i = 1, np
         pfit = pabs(i) - bsav(i)
         write(6,1013) elist(i),pabs(i),pfit, pbck(i),pres(i),bsav(i)
 15   continue
c
 1013 format ( 1x, f11.6, f16.6, 4f14.6 )
c
      rms = sqrt(phi/np)
c
      write(6,1014) rms
 1014 format ( '0rms error of fit = ', g15.3 )
c
      go to 1
c
 500  continue
      write(6,1500) ichan(chan), ityp(type), rsn, initial, final, np
 1500 format ( ' too few entries for ',/, 1x,a8,1x,a5,1x,a10,5x,
     >         3hi =, i5, 5h, f =, i5, " np = ", i2 )
      go to 1
c
 100  continue
      end







