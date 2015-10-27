*deck @(#)m7004.f	2.1  10/10/91
      program m7004
      implicit integer(a-z)
      common ia(1)
      real*8 z(1), charge, fpkey, rmin, rmax, rdel, energy, rswtch
      character *4096 ops
      character *16 cpass, fptoc
      character *3 itoc, ans
      character *800 card
      character *128 filci, filnm
      logical logkey, drven, tstspl, nospln, noireg
      equivalence (ia(1),z(1))
      common /memory/ ioff
      common/io/inp,iout
c
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "atomic ci filename" from rwf',-1,
     1             0,0,filci)
      call iosys ('open atomci as unknown',0,0,0,filci)
      level=intkey(ops,'print=m7004=level',0,' ')
      drven=logkey(ops,'m7004=driven-equation',.false.,' ')
      tstspl=logkey(ops,'m7004=test-spline',.false.,' ')
      nospln=logkey(ops,'m7004=no-spline',.false.,' ')
      noireg=logkey(ops,'m7004=no-irregular-function',.false.,' ')
      call posinp('$coul',cpass)
      call cardin(card)
      charge=fpkey(card,'charge',-1.d0,' ')
      rmin=fpkey(card,'minimum-r-value',1.d-20,' ')          
      rmax=fpkey(card,'maximum-r-value',10.d+00,' ')          
      rswtch=fpkey(card,'switching-r-value',.5d0,' ')
      rdel=fpkey(card,'integration-step-size',.01d0,' ')
      nener=intkey(card,'number-of-energies',1,' ')
      lmax=intkey(card,'maximum-angular-momentum',10,' ')
      ntrms=intkey(card,'maximum-number-terms-in-series',50,' ')
      toskp=intkey(card,'number-to-skip',4,' ')
c**********************************************************************c
c               set up parameters for integration                      c
c                            and                                       c
c                        interpolation                                 c
c               the integration uses a fine mesh                       c
c               which is then interpolated on a coarser grid           c     
c**********************************************************************c
      nint=(rmax-rmin)/rdel
      nint=nint+2
      nspln=0
      do 10 pt=1,nint,toskp
         nspln=nspln+1
   10 continue     
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
         ind=wpadti(scr+2*nint)
         c=iadtwp(ind)
         ftmp=c
         xtmp=ftmp
         words=ind
         if (.not.nospln) then 
             c=iadtwp(ind+nint)
             ftmp=c+nint
             xtmp=ftmp+nspln
             words=xtmp+nspln
         endif
         if (i.eq.1) then
             need=wpadti(words)
             call getscm(need,z,maxcor,'m7004',0)
         endif
   20 continue     
      write(iout,*) 'm7004:coulomb functions'
      write(iout,*) '                           print level = ',level
c
c
      write(iout,1) rmin,rmax,rdel,nint,rswtch
      if (.not.nospln) then
          write(iout,2) nspln, toskp
      endif
      if (noireg) then
          write(iout,3)
      endif
      do 30 ene=1,nener
         call posinp('$energy-'//itoc(ene),cpass)
         call cardin(card)         
         energy=fpkey(card,'energy',1.d0,' ')
         nol=intkey(card,'number-l-values',1,' ')
         call intarr(card,'l-values',ia(lval),nol,' ')
         if(.not.nospln) then
            filnm='"spline arrays-'//fptoc(energy)//'"'
            filsiz=nspln*nol*4
            call iosys ('does '//filnm//' exist on atomci',0,0,0,ans)
            if (ans.eq.'no') then
                call iosys ('create real '//filnm//' on atomci',
     1                       filsiz,0, 0,' ')
            endif
         endif
         call gencou(z(a),z(b),z(ak),z(bk),z(a0),z(b0),z(c0),z(d0),
     1               z(e0),z(psi0),z(psi1),z(dpsi0),z(dpsi1),charge,
     2               energy,z(v),z(x),z(xinv),z(c),z(scr),z(ftmp),
     3               z(xtmp),rmin,rmax,rdel,rswtch,ia(lval),nint,nol,
     4               toskp,nspln,filnm,nospln,noireg,ntrms,level)
         if (.not.nospln) then
             if (tstspl) then
                 call iosys ('rewind '//filnm//' on atomci '//
     1                       'read-and-write',0,0,0,' ')
                 call spltst(z(psi0),z(dpsi0),z(psi1),z(dpsi1),
     1                       z(space),z(x),z(ftmp),z(xtmp),z(c),
     2                       ia(ind),nint,nspln,filnm,nol)
             endif
         endif   
   30 continue     
      call iosys ('rewind all on atomci read-and-write',0,0,0,' ')
c
      call chainx(0)
c
    1 format(/,5x,'integration parameters:',/,20x,'r min =     ',e15.8,
     1         1x,'r max = ',e15.8,/,20x,'step size = ',e15.8,1x,
     2            'number steps = ',1x,i5,/20x,'series to integration '
     3            'switching point',1x,e15.8)           

    2 format(/,5x,'spline information:',/,20x,
     1            'number spline points = ',1x,i2,1x,
     2            'number to skip = ',1x,i2)
    3 format(/,20x,'no irregular solution computed')
      stop
      end

