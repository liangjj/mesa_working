      program resfit ( input, output, tape6, tape7, tape9, tape88 )
c
      integer chan,type,initial,final,rsn
      integer name(3), ityp(2), ichan(3)
      logical fit, debug, skip7, read88
      complex sfin(20)
      common / doyle / dum1, nit, dum2(22)
      common / kar   / dum4(150), bl(30),bu(30),v(30), dum5(60),
     >                 nuflag(30), nbd(30)
      common / kchg  / gsav(300,15), bsav(300), nersdm
      common / kvar  / npar, nprob, ners, dum3(8)
      common / pdata / np, npbck, elist(256), plist(256), sreal(256),
     >                 simag(256), phase(256), pabs(256), pres(256),
     >                 pbck(256), edup(256), ilist(256), del(256)
c
      namelist / inp / emin,emax, chan,type,initial,final, fit,debug,
     >                 nit,npbck,escal,eref,rsn,skip7,elist,del,name,
     >                 np, xnpi, read88
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
      data ityp / 8hinelastc, 8hreactive /
      data ichan / 5halpha, 5hbeta , 5hgamma /
c      data ebar / 50*0e0 /
c
      call filerep
c
      emin = -999e0
      emax = -999e0
      pi = acos(-1e0)
      chan = 1
      type = 2
      initial = 0
      final = -1
      fit = .true.
      debug = .false.
      de = 0e0
      ne = -1
      escal = 1e0
      eref = 0e0
      nbar = 0
      rsn = 10h
      npar = 3
      npbck = 2
      skip7 = .false.
      xnpi = 2.
      np = 0
      name(1) = 1h
      name(2) = 1h
      name(3) = 1h
      nit = 10
      rsn = 1h
      read88 = .false.
c
      nread = 0
 1    continue
      read inp
      if ( eof(5hinput) ) 100, 2
 2    continue
c
      do 3 i = 1, npar
      read *, nuflag(i), v(i), nbd(i),bl(i),bu(i)
 3    continue
      ntpar = npar + npbck + 1
c
      if ( skip7 ) go to 50
c
c     scan through tape7 file for transition desired
c
      rewind 7
      nread = nread + 1
      nc = 0
      ip = 0
c
 4    continue
c
      read ( 7, 1001 ) it,ic,name,etot,is,psum,nfs,irsn
 1001 format ( 2x,a8,1x,a5,5x,3a2,4x,e20.12,i10,e20.12,i10,a9 )
      if ( eof(7) ) 6, 5
 5    nc = nc + 1
      read ( 7, 1002 ) (sfin(i),i=1,nfs)
 1002 format ( 2x, 10e11.3 )
c
c     check for proper input card
c
      if ( irsn .ne. rsn           ) go to 4
      if (  it  .ne. ityp(type)    ) go to 4
      if (  ic  .ne. ichan(chan)   ) go to 4
      if (  is  .ne. initial       ) go to 4
      if ( nfs  .lt. 1+final       ) go to 4
      if ( etot .lt. emin          ) go to 4
      if ( etot .gt. emax          ) go to 4
c
      ifs = 1 + max0(final,0)
      pfin = sfin(ifs)*conjg(sfin(ifs))
      ip = ip + 1
      elist(ip) = escal*(etot - eref)
      sreal(ip) = real(sfin(ifs))
      simag(ip) = aimag(sfin(ifs))
      plist(ip) = pfin
      if ( final .lt. 0 ) plist(ip) = psum
c
      go to 4
c
c     sort entries by energy
c
 6    continue
      np = ip
      if ( np .lt. ntpar ) go to 500
c
      call fsort1 ( np, 3, elist, ilist, plist, sreal, simag )
c
c     discard duplicates
c
      i = 1
      idup = 0
    7 continue
      ip1 = i + 1
      if ( ip1 .gt. np ) go to 9
      if ( abs(elist(ip1)-elist(i)) .lt. 1e-7 ) go to 8
      i = ip1
      go to 7
 8    continue
      idup = idup + 1
      edup(idup) = elist(i)
      call condens ( elist, plist, sreal, simag, i, np )
      go to 7
    9 continue
c
      do 10 i = 1, np
      phase(i) = atan2(simag(i), sreal(i))/pi
 10   continue
      go to 51
c
c     input of phases given directly
c
 50   continue
      if ( read88 ) read(88,*) (ipt,elist(i),del(i),i=1,np)
      idup = 0
      do 52 i = 1, np
      sreal(i) = 0e0
      simag(i) = 0e0
      plist(i) = 0e0
      phase(i) = del(i)/pi
 52   continue
c
 51   continue
      call phasit ( elist, phase, np, xnpi, 0.,-1., pres, pabs )
      do 11 i = 1, np
      pabs(i) = pabs(i)*pi
 11   continue
c
      print 1300, name,ityp(type),ichan(chan),initial,final,rsn
 1300 format ( "1sorted list of s-matrix elements for the transition",/,
     >         10x, 3a2,1x,a8,1x,a5, " from v=",i2,"    to v="i2,
     >         ",   rsn=",a10, /,
     >          14x, "energy     smat**2     s(real)     s(imag)",
     >          "       phase/pi    abs phase" )
      print 1301, (elist(i),plist(i),sreal(i),simag(i),phase(i),pabs(i),
     >            i=1,np)
 1301 format ( 5x, 0pf15.7, 1p3e12.3, 0pf15.4, f13.4 )
      if ( idup .gt. 0 ) print 1302, idup, (edup(i),i=1,idup)
 1302 format ( * a total of *, i3,* duplicate energies found -- *,/,
     >         (1x,10f13.6) )
c
      print 1010, npbck, nit,
     >            (nuflag(i), v(i), nbd(i), bl(i), bu(i), i=1, 3)
 1010 format ( "1   non-linear least squares resonance fitting program",
     >         /,"0resonant part of phase shift fitted to 3 nonlinear"
     >           " parameters",
     >         /,"      pres(e)= p*(.5*pi + atan(2*(e-eres)/gres))",
     >         /,"          parameters are -- p,eres,gres",
     >         /,"0background part of phase shift fitted to e**",i2,
     >         /,"0maximum number of iterations = ", i3,/,
     >         /,"0",12x,"vary flag    initial guess    bound flag",
     >           "    lower bound    upper bound",
     >         /,"0   -p-      ", i6, g20.6, i10, g19.6, g15.6,
     >         /,"   -eres-    ", i6, g20.6, i10, g19.6, g15.6,
     >         /,"   -gres-    ", i6, g20.6, i10, g19.6, g15.6 )
c
      ners = np
      call opac
      call fcn ( v, bsav, phi )
c
      print 1011, (v(i),i=1,ntpar)
 1011 format ( "1final fit of data follows --", /,
     >        "0resonant amplitude      p    = ", g15.5, /,
     >        "0resonance energy     eres    = ", g15.5, /,
     >        "0resonance width      gres    = ", g15.5, /,
     >        "0linear background parameters = ", 5g15.5,/,(30x,5g15.5))
      print 1012
 1012 format ( //,"0     energy    actual phase    fit phase    ",
     >         "background     resonance         error" )
      do 15 i = 1, np
      pfit = pabs(i) - bsav(i)
      print 1013, elist(i),pabs(i),pfit, pbck(i),pres(i),bsav(i)
 15   continue
c
 1013 format ( 1x, f11.6, f16.6, 4f14.6 )
c
      rms = sqrt(phi/np)
c
      print 1014, rms
 1014 format ( "0rms error of fit = ", g15.3 )
c
      go to 1
c
 500  continue
      print 1500, ichan(chan), ityp(type), rsn, initial, final, np
 1500 format ( * too few entries for *,/, 1x,a8,1x,a5,1x,a10,5x,
     >         3hi =, i5, 5h, f =, i5, " np = ", i2 )
      go to 1
c
 100  continue
      end
