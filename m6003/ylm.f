*deck @(#)ylm.f	1.1 9/8/91
c***begin prologue     m6003
c***date written       881006   (yymmdd)
c***revision date      890418   (yymmdd)
c***keywords           m6003, link 6003, spherical harmonics
c***author             rescigno, t. n.(llnl)
c***source             m6003
c***purpose            spherical harmonic generation
c***description        calculates spherical harmonics on grid
c***                   for integral calculation. fully vectorized
c*** 
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m6003
      program ylm 
      implicit integer (a-z)
      character *3 ans
      character *8 cpass, chrkey
      character *128 filgrd, filsph
      character *1600 card
      character *4096 ops
      character *13 grdnam
      logical grdtyp, logkey, prnt
      real *8 z
      dimension z(1)
      common a(1)
      common /io/ inp, iout
      common /memory / ioff
      equivalence (z,a)
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      prnt=logkey(ops,'print=m6003=ylm',.false.,' ')
      write (iout,100)
      call posinp ('$ylm',cpass)
      call cardin(card)
      grdtyp=logkey(card,'untransformed-grid',.false.,' ')
      grdnam='"trns grid"'
      if (grdtyp) then
          grdnam='"untrns grid"'
      endif
      lmax=intkey(card,'maximum-l-value',30,' ')
      mumax=intkey(card,'maximum-m-value',5,' ')
      write (iout,300) lmax, mumax
      call iosys ('read character "grid filename" from rwf',-1,0,0,
     1             filgrd)
      call iosys ('read character "spherical harmonic filename" '//
     1            'from rwf',-1,0,0,filsph)
      call iosys ('open grid as old',0,0,0,filgrd)
      call iosys ('read integer "no. grid pts" from grid',1,ngrid,0,
     1            ' ')
      call iosys ('read integer "point buffer" from grid',1,
     1            pntbuf,0,' ')
      write (iout,400) ngrid, pntbuf
      call iosys ('open ylms as new on ssd',262144,0,0,filsph)
      call iosys ('write character "grid type" to ylms',0,0,0,grdnam)
      call iosys ('write integer "no. grid pts" to ylms',1,ngrid,0,
     1            ' ')
      call iosys ('write integer "point buffer" to ylms',1,pntbuf,0,
     1            ' ')
      call iosys ('write integer "max l in ylm" to ylms',1,lmax,0,' ')
      call iosys ('write integer "max m in ylm" to ylms',1,mumax,0,' ')
      call iosys ('create real ylm on ylms',-1,0,0,' ')
c----------------------------------------------------------------------c
c              calculate number of passes needed                       c
c----------------------------------------------------------------------c
      npass=ngrid/pntbuf
      nwleft=ngrid-pntbuf*npass
      if (nwleft.ne.0) then
          npass=npass+1
      else
          nwleft=pntbuf
      endif
      lmmax=lmax+mumax
c----------------------------------------------------------------------c
c                     memory allocation                                c
c----------------------------------------------------------------------c
      call iosys ('read integer maxsiz from rwf',1,canget,0,' ')
      maxfac=lmmax+10
      words=2*maxfac+11*pntbuf+3*pntbuf*(lmax+1)
      if (canget.lt.words) then
          call lnkerr('not enough memory:decrease point buffer')
      endif
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(words,z,ngot,'m6003',0)      
      dfct=ioff
      ddfct=dfct+maxfac
      grid=ddfct+maxfac
      x=grid+4*pntbuf
      cphi=x+pntbuf
      sphi=cphi+pntbuf
      plm=sphi+pntbuf
      cmp=plm+pntbuf*(lmax+1)
      cmphi=cmp+pntbuf+pntbuf
      smphi=cmphi+pntbuf
      ylmp=smphi+pntbuf
      ylmm=ylmp+pntbuf*(lmax+1)
      words=ylmm+pntbuf*(lmax+1)
 
      call fact(z(dfct),z(ddfct),maxfac)
c----------------------------------------------------------------------c
c                    loop over points                                  c
c----------------------------------------------------------------------c
      do 10 ipass=1,npass
         nwread=pntbuf
         if (ipass.eq.npass) then
             nwread=nwleft
         endif
c----------------------------------------------------------------------c
c                    read buffer of points                             c
c----------------------------------------------------------------------c        
         call iosys ('read real '//grdnam//' from grid without '//
     1               'rewinding',4*nwread,z(grid),0,' ')
c----------------------------------------------------------------------c
c                 calculate some necessary functions                   c
c                 used over and over again in plm routine              c
c----------------------------------------------------------------------c 
         call miscfn(z(grid),z(x),z(cphi),z(sphi),nwread)
c----------------------------------------------------------------------c
c             calculate the legendre functions                         c
c----------------------------------------------------------------------c
         do 20 mu=0,mumax
            call legend(z(plm),z(x),z(dfct),z(ddfct),nwread,lmax,mu,
     1                  maxfac)
            call sph(z(plm),z(cphi),z(sphi),z(cmp),z(cmphi),
     1               z(smphi),z(ylmp),z(ylmm),z(dfct),nwread,mu,lmax,
     2               maxfac,ipass)
            if (prnt) then
                call prnty(z(plm),z(x),nwread,mu,lmax)
            endif
   20    continue
   10 continue
c----------------------------------------------------------------------c
c              close random file                                       c
c----------------------------------------------------------------------c
      call iosys ('endfile ylm on ylms',0,0,0,' ')
      call iosys ('rewind all on ylms read-and-write',0,0,0,' ')
      call iosys ('rewind all on grid read-and-write',0,0,0,' ')
      call iosys('close ylms',0,0,0,' ')
      call iosys ('close grid',0,0,0,' ')
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)
      stop
  100 format (//,20x,'***** m6003: calculate spherical harmonics *****')
  300 format(/,5x,'maximum l value',1x,i3,5x,'maximum m value',1x,i3)
  400 format(/,5x,'no. points',1x,i8,5x,'point buffer',1x,i8)
      end
