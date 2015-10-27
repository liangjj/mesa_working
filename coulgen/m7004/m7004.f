*deck @(#)m7004.f	2.1  10/10/91
      program m7004
      implicit integer(a-z)
      parameter ( nen=100 )
      common ia(1)
      real*8 z(1), charge, fpkey, rmin, rmax, rdel, energy, rswtch
      real*8 alpha
      dimension energy(nen)
      character *4096 ops
      character *16 cpass, fptoc
      character *3 itoc, ans
      character *800 card
      character *128 filci, filnm
      character *80 title
      logical logkey, drven, tstspl, nospln, noireg, prgrd, prser, prfun
      equivalence (ia(1),z(1))
      common /memory/ ioff
      common/io/inp,iout
c
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "atomic ci filename" from rwf',-1,
     1             0,0,filci)
      call iosys ('open atomci as new',0,0,0,filci)
      drven=logkey(ops,'m7004=driven-equation',.false.,' ')
      tstspl=logkey(ops,'m7004=test-spline',.false.,' ')
      nospln=logkey(ops,'m7004=no-spline',.false.,' ')
      noireg=logkey(ops,'m7004=no-irregular-function',.false.,' ')
      prgrd=logkey(ops,'print=m7004=grid',.false.,' ')
      prser=logkey(ops,'print=m7004=series',.false.,' ')
      prfun=logkey(ops,'print=m7004=functions',.false.,' ')
      call posinp('$coul',cpass)
      call cardin(card)
      charge=fpkey(card,'charge',-1.d0,' ')
      rmin=fpkey(card,'minimum-r-value',1.d-20,' ')          
      rmax=fpkey(card,'maximum-r-value',10.d+00,' ')          
      alpha=fpkey(card,'alpha-value',1.d0,' ')          
      rswtch=fpkey(card,'switching-r-value',.5d0,' ')
      rdel=fpkey(card,'integration-step-size',.01d0,' ')
      order=intkey(card,'order-of-spline-fit',4,' ')
      nener=intkey(card,'number-of-energies',1,' ')
      lmax=intkey(card,'maximum-angular-momentum',10,' ')
      ntrms=intkey(card,'maximum-number-terms-in-series',50,' ')
      toskp=intkey(card,'number-to-skip',4,' ')
      call iosys('write real "smallest r" to atomci',1,rmin,0,' ')
      call iosys('write real "largest r" to atomci',1,rmax,0,' ')
      call iosys('write real "step size" to atomci',1,rdel,0,' ')
      call iosys('write real "switching r" to atomci',1,rswtch,0,' ')
      call iosys('write integer "largest l" to atomci',1,lmax,0,' ')
      call iosys('write integer "number of terms in series" to atomci',
     1            1,ntrms,0,' ')
c**********************************************************************c
c               set up parameters for integration                      c
c                            and                                       c
c                        interpolation                                 c
c               the integration uses a fine mesh                       c
c               which is then interpolated on a coarser grid           c     
c**********************************************************************c
c**********************************************************************c
c    the next bit of fancy footwork is to make sure in a simpson       c
c    rule integration we end on the first and last points in the       c
c                           integration                                c
c**********************************************************************c
      nint=(rmax-rmin)/rdel+1
      panel=(nint-4)/3 +1
      nint=3*panel+4
c**********************************************************************c
c              find largest x for series expansion                     c
c**********************************************************************c
      call fninpt(rswtch,rmin,rdel,nint,ptbeg)
      write(iout,5) ptbeg
c**********************************************************************c
c                spline fit is only done for integration points        c
c**********************************************************************c
      nspln=0
      if (.not.nospln) then
          do 10 pt=ptbeg,nint,toskp
             nspln=nspln+1
   10     continue     
          call iosys('write integer "number of spline points" to '//
     1               'atomci',1,nspln,0,' ')
          nbreak=nspln-order+1
          call iosys('write integer "number of break points" to atomci',
     1                1,nbreak,0,' ')
      endif
      ioff=1
      do 20 i=1,2
         lval=wpadti(ioff)
         a=iadtwp(lval+lmax+1)
         b=a+ntrms+1
         ak=b+ntrms+1
         bk=ak+ntrms+1
         a0=bk+ntrms+1
         b0=a0+ntrms+1
         c0=b0+ntrms+1
         d0=c0+ntrms+1
         e0=d0+ntrms+1
         xinv=e0+ntrms+1
         v=xinv+nint
         x=v+nint
         space=x
         psi0=x+nint
         dpsi0=psi0+nint
         psi1=dpsi0+nint
         dpsi1=psi1+nint
         scr=dpsi1+nint
         ind=wpadti(scr+nint)
         break=iadtwp(ind)
         c=break
         sc=c
         ftmp=sc
         xtmp=ftmp
         rhs=xtmp
         sdrv=rhs
         words=scr+nint
         if (.not.nospln) then 
             ind=wpadti(scr+nint)
             break=iadtwp(ind+nint)
             c=break+nbreak+1
             sc=c+order*nbreak
             ftmp=sc+(4*nspln+1)*order*2
             xtmp=ftmp+nspln
             words=xtmp+nspln
         endif
         if (drven) then
             rhs=words
             sdrv=rhs+nint
             words=sdrv+nint
         endif 
         if (i.eq.1) then
             need=wpadti(words)
             call getscm(need,z,maxcor,'m7004',0)
         endif
   20 continue     
      write(iout,*) 'm7004:coulomb functions'
c
c
c**********************************************************************c
c                    make integration and spline grid                  c
c**********************************************************************c
      write(iout,1) rmin,rmax,rdel,nint,rswtch
      if (.not.nospln) then
          write(iout,2) order,toskp,nspln
      endif
      if (noireg) then
          write(iout,3)
      endif
      if (drven) then
          write(iout,*)
          write(iout,*) '                   solve driven equation'
      endif
      call mkgrd(z(x),z(xinv),z(xtmp),rmin,rmax,rdel,nint,ptbeg,last,
     1           nospln,toskp)
      if (prgrd) then
          title='integration grid'
          write(iout,*) title
          write(iout,4) (z(i),i=x,x+nint-1)
      endif
      if (drven) then
          call mkexp(z(x),z(rhs),alpha,nint)
      endif
c**********************************************************************c
c           calculate the spline information only dependent on the     c
c           grid used not on the fitted function.                      c
c**********************************************************************c
      if (.not.nospln) then
          call prespl(nspln,z(xtmp),order,z(sc))
          call splmat(nspln,z(xtmp),order,z(sc))
      endif
      do 30 ene=1,nener
         call posinp('$energy-'//itoc(ene),cpass)
         call cardin(card)         
         energy(ene)=fpkey(card,'energy',1.d0,' ')
         nol=intkey(card,'number-l-values',1,' ')
         call intarr(card,'l-values',ia(lval),nol,' ')
         call iosys('write integer "number of l values energy-'
     1               //fptoc(energy(ene))//'" to atomci',1,nol,0,' ')
         call iosys('write integer "l values energy-'//
     1              fptoc(energy(ene))//'" to atomci',nol,
     2              ia(lval),0,' ')
         call gencou(z(a),z(b),z(ak),z(bk),z(a0),z(b0),z(c0),z(d0),
     1               z(e0),z(psi0),z(psi1),z(dpsi0),z(dpsi1),charge,
     2               energy(ene),z(v),z(x),z(xinv),z(break),z(c),z(sc),
     3               z(scr),z(ftmp),z(xtmp),z(rhs),z(sdrv),rmin,rmax,
     4               rdel,rswtch,ptbeg,last,ia(lval),nint,nol,order,
     5               nbreak,toskp,nspln,filnm,nospln,noireg,ntrms,
     6               drven,prser,prfun)
   30 continue
      if (.not.nospln) then
          call iosys ('write real breakpoints to atomci',nbreak+1,
     1                 z(break),0,' ')     
          if (tstspl) then
              do 40 ene=1,nener
                 call extcou(z(psi0),z(dpsi0),z(psi1),z(dpsi1),z(a),
     1                       z(b),z(a0),z(b0),z(c0),z(space),z(x),
     2                       z(xinv),z(break),z(c),ia(ind),ia(lval),
     3                       energy(ene),nint,ntrms,order,nbreak,
     4                       toskp,filnm,nol)
   40         continue
          endif
      endif   
      call iosys ('rewind all on atomci read-and-write',0,0,0,' ')
c
      call chainx(0)
c
    1 format(/,5x,'integration parameters:',/,20x,'r min =     ',e15.8,
     1         1x,'r max = ',e15.8,/,20x,'step size = ',e15.8,1x,
     2            'number steps = ',1x,i5,/20x,'series to integration '
     3            'switching point',1x,e15.8)           

    2 format(/,5x,'spline parameters:',/,20x,'spline order = ',1x,i2,1x,
     1            'number to skip = ',1x,i2,1x,
     2            'number spline points = ',1x,i5)
    3 format(/,20x,'no irregular solution computed')
    4 format(5x,5e15.8)   
    5 format(/,5x,'integrate differential equation from point',1x,i5,
     1            1x,'onward')
      stop
      end

