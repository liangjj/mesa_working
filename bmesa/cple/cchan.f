*deck cchan 
c***begin prologue     m6235
c***date written       930627   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6235, link 6235, bessel
c***author             schneider, b. i.(nsf)
c***source             m6235
c***purpose            coupled channel r-matrix code using ricatti-bessel 
c***                   functions as translational basis.
c***description        a mixed basis of ricatti-bessel functions having      
c***                   zero function value and zero derivative value are
c***                   used to solve coupled channel equations. 
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m6235
      program cchan
c
      implicit integer (a-z)
      parameter ( maxrt=200, maxl=10, ptmax=1000, chnmax=20 )
      parameter ( noene=200 )
      dimension nobf(chnmax), nobfz(chnmax), nobfzd(chnmax)
c      dimension nobfex(chnmax)
      dimension lval(chnmax), ec(chnmax), wt(ptmax), point(ptmax)
      dimension root(maxrt,chnmax), droot(maxrt,chnmax), energy(noene)
      dimension vcij(chnmax*(chnmax+1)/2)
      character*4096 ops
      character*10 cpass
      character*1600 card
      character*128 fillam
      character*3 itoc, chr1, chr2
      character*32 chrkey, pottyp
      character*80 title
      logical posinp, logkey, prntwv, prntgr, prntir, dopse
      logical prntov, prntke, prntpe, prnth, prntld, prntkm, prntep
      logical elist, preigs, scat
      common z(1)
      dimension ia(1)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
      real*8 z, rbox, fpkey, root, droot, wt, point, rst, rfn, rfs
      real*8 ec, tmp, wtsum, e0, deltae, energy, psesum, vcij, range
      real*8 rfinal
c
      call drum
      write(iout,*)
      write(iout,*) '                         m6235:coupled channel'//
     1                                        ' r-matrix link'
      call iosys ('read character options from rwf',-1,0,0,ops)
      prntwv=logkey(ops,'print=m6235=wavefunction',.false.,' ')
      prntgr=logkey(ops,'print=m6235=grid',.false.,' ')
      prntov=logkey(ops,'print=m6235=overlap-matrix',.false.,' ')
      prntke=logkey(ops,'print=m6235=kinetic-energy-matrix',.false.,' ')
      prntpe=logkey(ops,'print=m6235=potential-energy-matrix',
     1              .false.,' ')
      preigs=logkey(ops,'print=m6235=hamiltonian-eigenvalues',
     1                   .false.,' ')
      prnth=logkey(ops,'print=m6235=hamiltonian-matrix',.false.,' ')
      prntld=logkey(ops,'print=m6235=r-matrix',.false.,' ')
      prntkm=logkey(ops,'print=m6235=k-matrix',.false.,' ')
      prntep=logkey(ops,'print=m6235=eigenphase-vectors',.false.,' ')
      prntir=logkey(ops,'print=m6235=intermediate-results',
     1              .false.,' ')
      pottyp=chrkey(ops,'potential-type','square-well',' ')
      scat=logkey(ops,'m6235=scattering=on',.false.,' ')
      call iosys ('read character "linear algebraic filename" from rwf',
     1            -1,0,0,fillam)
      call iosys ('open lamdat as old',0,0,0,fillam)
      if (posinp('$grid',cpass) ) then
          call cardin(card)
      endif          
      call iosys ('write character "potential type" to lamdat',0,0,
     1             0,pottyp)
      rbox=fpkey(card,'r-matrix-radius',10.d0,' ')
      rfinal=fpkey(card,'final-r',rbox,' ')
      call iosys('write real "r-matrix radius" to lamdat',1,rbox,0,' ')
      call iosys('write real "final radius" to lamdat',1,rfinal,0,' ')
      range=fpkey(card,'potential-range',rbox,' ')
      call iosys('write real "potential range" to lamdat',1,
     1            range,0,' ')
      nreg=intkey(card,'number-of-integration-regions',1,' ')
      ptcnt=0
      wtsum=0.d0
      do 10 ireg=1,nreg
         npn=intkey(card,'number-of-points-region-'//itoc(ireg),4,' ')
         rst=fpkey(card,'starting-r-region-'//itoc(ireg),0.d0,' ')
         rfn=fpkey(card,'final-r-region-'//itoc(ireg),rbox,' ')
         write(iout,1) ireg, npn, rst,rfn
         rfs=rfn-rst
         if (prntgr) then
             write(iout,2)
         endif
         do 20 pt=1,npn
            ptcnt=ptcnt+1
            call lgndrx (npn,pt,wt(ptcnt),point(ptcnt))
            wt(ptcnt)=wt(ptcnt)*rfs
            wtsum=wtsum+wt(ptcnt)
            point(ptcnt)=point(ptcnt)*rfs+rst
            if (prntgr) then
                write(iout,3) ptcnt, point(ptcnt), wt(ptcnt)
            endif                
   20    continue
   10 continue
      ptcnt=ptcnt+1
      wt(ptcnt)=0.d0
      point(ptcnt)=rbox
      write(iout,7) wtsum
      if ( posinp('$channels',cpass) ) then
           call cardin(card)
      endif
      nc=intkey(card,'number-of-channels',1,' ')
      call iosys ('write integer "number of channels" to lamdat',
     1             1,nc,0,' ')
      call intarr(card,'channel-l-values',lval,nc,' ')
      call iosys ('write integer "channel l values" to lamdat',nc,
     1             lval,0,' ')              
      call fparr(card,'channel-energies',ec,nc,' ')
      call iosys ('write real "channel energies" to lamdat',nc,
     1             ec,0,' ')      
      call intarr(card,'number-of-basis-functions-per-channel',nobf,
     1            nc,' ')
      call intarr(card,'number-of-zero-value-basis-functions-per-'//
     1            'channel',nobfz,nc,' ')
      call intarr(card,'number-of-zero-derivative-basis-functions-'//
     1            'per-channel',nobfzd,nc,' ')
      ij=0
      do 30 i=1,nc
         call pakstr(itoc(i),ii)
         do 40 j=1,i
            ij=ij+1
            call pakstr(itoc(j),jj)
            chr1=itoc(i)
            chr2=itoc(j)
            vcij(ij)=fpkey(card,'v'//chr1(1:ii)//chr2(1:jj),-1.d0,' ')
 40      continue
 30   continue
      ncc=0
      nco=0
      mxchn=0
      do 50 i=1,nc
c         nobfex(i)=nobf(i)-nobfz(i)-nobfzd(i)
         nobfzd(i)=nobf(i)-nobfz(i)
         mxchn=max(mxchn,nobf(i))
c         ncc=ncc+nobfz(i)+nobfex(i)
         ncc=ncc+nobfz(i)
         nco=nco+nobfzd(i)
         call iosys ('read integer "number of j roots '//
     1                    'for l='//itoc(lval(i))//'" from lamdat',
     2                     1,nozrt,0,' ')
         if (nozrt.lt.nobfz(i)) then
             call lnkerr('error in number of zero function roots')
         else
             call iosys ('read real "j roots for '//
     1                        'l='//itoc(lval(i))//'" from lamdat',
     2                         nobfz(i),root(1,i),0,' ')
             do 60 iroot=1,nobfz(i)
                root(iroot,i)=root(iroot,i)/rbox
   60        continue                  
         endif
         call iosys ('read integer "number of jp roots '//
     1               'for l='//itoc(lval(i))//'" from lamdat',1,
     2                nozdrt,0,' ')
         if (nozdrt.lt.nobfzd(i)) then
             call lnkerr('error in number of zero derivative '//
     1                   'roots')
         else
             call iosys ('read real "jp roots for '//
     1                   'l='//itoc(lval(i))//'" from lamdat',
     2                    nobfzd(i),droot(1,i),0,' ')
             do 70 iroot=1,nobfzd(i)
                droot(iroot,i)=droot(iroot,i)/rbox
   70        continue         
         endif
   50 continue                
      matsiz=ncc+nco
      write(iout,4) nc, rbox
      write(iout,5)
      lmax=0
      do 80 i=1,nc
         lmax=max(lmax,lval(i))
         ii=i
c         write(iout,6) i, lval(i), nobfz(i), nobfzd(i), nobfex(i), ec(i)
         write(iout,6) i, lval(i), nobfz(i), nobfzd(i), ec(i)
 80   continue
      call iosys ('write integer "maximum l value" to lamdat',1,
     1             lmax,0,' ')     
      tmp=lmax*30.d0
      ltop=lmax+sqrt(tmp)
      ltop=max(ltop,lmax)
      lp=ltop+1     
      ioff=1
      do 90 i=1,2
         f=ioff
         df=f+ptcnt*matsiz
         ddf=df+ptcnt*matsiz
         flst=ddf+ptcnt*matsiz
         dflst=flst+matsiz 
         x=dflst+matsiz 
         j=x+ptcnt
         jp=j+lp*ptcnt
         y=jp+lp*ptcnt
         yp=y+lp*ptcnt
         wron=yp+lp*ptcnt
         scr=wron+lp
         w1=scr+2*ptcnt
         w1=wpadti(w1)
         pot=x
         s=pot+ptcnt
         t=s+mxchn*mxchn
         v=t+mxchn*mxchn*nc
         eig=v+mxchn*mxchn*nc*nc
         dum=eig+matsiz
         scr1=dum+matsiz
         ham=scr1+max(2,matsiz,ptcnt*mxchn)
         gamic=ham+matsiz*matsiz
         rmat=gamic+nc*matsiz
         words=wpadti(rmat+nc*nc)
         if (scat) then
             jbox=rmat+nc*nc
             jpbox=jbox+lp
             ybox=jpbox+lp
             ypbox=ybox+lp
             wronbx=ypbox+lp
             y0=wronbx+lp
             dy0=y0+nc
             y1=dy0+nc
             dy1=y1+nc
             coef=dy1+nc
             kmat=coef+2*nc*nc
             kmatc=coef
             tmtrx=kmat+nc*nc
             cross=tmtrx+2*nc*nc
             ipvt=wpadti(cross+nc*nc)
             w2=ipvt+nc
             words=max(w1,w2)
         endif             
         if (i.eq.1) then
             call iosys ('read integer maxsiz from rwf',1,
     1                    canget,0,' ')
             if (words.gt.canget) then
                 call lnkerr('not enough memory. will quit')
             endif
             call iosys ('write integer maxsiz to rwf',1,
     1                    words,0,' ')
         else
             call getscm(words,z,ngot,'cple',0)
         endif
   90 continue
c      calculate the channel basis with zero value at r-matrix
c      boundary and any exponential basis functions if present.
c      these will be first in the list of r-matrix basis functions and
c      we call them the closed channel orbitals.
      f0=f
      df0=df
      ddf0=ddf
      if (ncc.ne.0) then
          do 100 i=1,nc
             call mkcfun(z(f0),z(df0),z(ddf0),point,root(1,i),z(x),
     1                   z(j),z(jp),z(y),z(yp),z(wron),z(scr),
     2                   lval(i),nobfz(i),ptcnt,prntwv)
             f0=f0+ptcnt*nobfz(i)
             df0=df0+ptcnt*nobfz(i)
             ddf0=ddf0+ptcnt*nobfz(i)
c             if (nobfex(i).ne.0) then
c                 call mkexp(z(f0),z(df0),z(ddf0),point,lval(i),
c     1                      nobfex(i),ptcnt,i,prntwv)
c                 f0=f0+ptcnt*nobfex(i)
c                 df0=df0+ptcnt*nobfex(i)
c                 ddf0=ddf0+ptcnt*nobfex(i)
c             endif
 100      continue   
      endif
c      calculate the channel basis with zero derivative at the
c      r-matrix boundary.  these will follow all of the basis functions
c      with zero value and are called the open channel orbitals.
      if (nco.ne.0) then
          do 200 i=1,nc
             call mkcfun(z(f0),z(df0),z(ddf0),point,droot(1,i),
     1                   z(x),z(j),z(jp),z(y),z(yp),z(wron),
     2                   z(scr),lval(i),nobfzd(i),ptcnt,prntwv)
             f0=f0+ptcnt*nobfzd(i)
             df0=df0+ptcnt*nobfzd(i)
             ddf0=ddf0+ptcnt*nobfzd(i)
 200      continue
      endif
      fi=f
      dfi=df
      ddfi=ddf
      fj=fi+ncc*ptcnt
      dfj=dfi+ncc*ptcnt
      ddfj=ddfi+ncc*ptcnt
      flsti=flst
      dflsti=dflst
      flstj=flsti+ncc
      dflstj=dflsti+ncc
      call prep(z(fi),z(dfi),z(ddfi),z(fj),z(dfj),z(ddfj),z(flsti),
     1          z(flstj),z(dflsti),z(dflstj),z(scr),wt,ptcnt,ncc,nco)
c      calculate the overlap and kinetic energy matrix elements.
c      they are diagonal in the channel index.
      t0=t
      newmat=0
      do 300 i=1,nc
         nbf=nobfz(i)+nobfzd(i)
c        calculate the overlap matrix of the primitive functions
         call smat(z(fi),z(fj),z(s),nobfz(i),nobfzd(i),ptcnt,
     1              nbf,prntir)
c        normalize the overlap
         call renrm(z(s),z(fi),z(dfi),z(ddfi),z(fj),z(dfj),z(ddfj),
     1              z(flsti),z(flstj),z(dflsti),z(dflstj),
     2              z(scr1),nobfz(i),nobfzd(i),ptcnt,nbf,prntir)
c        the closed channel space is an orthonormal space
c        orthogonalize the open channel orbitals to the closed chanel
c        orbitals and then transform the open channel functions to
c        an orthonormal basis.  in the end the entire set is orthonormal
         ni=nobfz(i)
         nj=nobfzd(i)
         nbf=ni+nj
         if (ni.ne.0.and.nj.ne.0) then 
             call toorth(z(fi),z(dfi),z(ddfi),z(fj),z(dfj),z(ddfj),
     1                   z(flsti),z(flstj),z(dflsti),z(dflstj),z(s),
     2                   z(eig),z(dum),z(scr1),nobfz(i),nobfzd(i),
     3                   ptcnt,nbf,nj,prntir)
             nbf=nobfz(i)+nj
         endif
         nobfz(i)=ni
         nobfzd(i)=nj
         newmat=newmat+ni+nj
         call smat(z(fi),z(fj),z(s),ni,nj,ptcnt,nbf,prntir)
         call tmat(z(fi),z(dfi),z(ddfi),z(fj),z(dfj),z(ddfj),z(flsti),
     1             z(flstj),z(dflsti),z(dflstj),z(t0),ni,nj,
     2             ptcnt,nbf,prntir)
         if (prntov) then
             title='overlap matrix for channel = '//itoc(i)
             call prntrm(title,z(s),nbf,nbf,nbf,nbf,iout)
         endif
         if (prntke) then
             title='kinetic energy matrix for channel = '//itoc(i)
             call prntrm(title,z(t0),nbf,nbf,nbf,nbf,iout)
         endif
         fi=fi+ptcnt*nobfz(i)
         dfi=dfi+ptcnt*nobfz(i)
         ddfi=ddfi+ptcnt*nobfz(i)
         fj=fj+ptcnt*nobfzd(i)
         dfj=dfj+ptcnt*nobfzd(i)
         ddfj=ddfj+ptcnt*nobfzd(i)
         flsti=flsti+nobfz(i)
         dflsti=dflsti+nobfz(i)
         flstj=flstj+nobfzd(i)
         dflstj=dflstj+nobfzd(i)
         t0=t0+nbf*nbf
  300 continue
      fci=f
      foi=fci+ncc*ptcnt
      vij=v
      ntri=nc*(nc+1)/2
      do 400 i=1,nc
         nbfi=nobfz(i)+nobfzd(i)
         fcj=f
         foj=fcj+ncc*ptcnt
         do 500 j=1,i
            nbfj=nobfz(j)+nobfzd(j)              
            call potntl(z(fci),z(foi),z(fcj),z(foj),vcij,z(pot),
     1                  z(vij),point,range,nobfz(i),nobfzd(i),nobfz(j),
     2                  nobfzd(j),ptcnt,nbfi,nbfj,pottyp,i,j,
     3                  ntri,prntir)
            if (prntpe) then
                call pakstr(itoc(i),ii)
                call pakstr(itoc(j),jj)
                chr1=itoc(i)
                chr2=itoc(j)
                title='potential energy matrix for channels = ('
     1                //chr1(1:ii)//','//chr2(1:jj)//')'
                call prntrm(title,z(vij),nbfi,nbfj,nbfi,nbfj,iout)
            endif                       
            fcj=fcj+nobfz(j)*ptcnt
            foj=foj+nobfzd(j)*ptcnt
            vij=vij+nbfi*nbfj
  500    continue
         fci=fci+nobfz(i)*ptcnt
         foi=foi+nobfzd(i)*ptcnt
  400 continue
      call filham(z(ham),z(t),z(v),ec,nobfz,nobfzd,nc,ncc,newmat)
      if (prnth) then
          title='final hamiltonian matrix'
          call prntrm(title,z(ham),newmat,newmat,newmat,newmat,iout)
      endif           
      call diag(z(ham),z(eig),z(dum),newmat,preigs)
      call surpsi(z(ham),z(flst),z(gamic),nobfz,nobfzd,nc,ncc,
     1            newmat,prntir)
      dopse=.false.
      if (scat ) then
          if (posinp('$energy',cpass)) then
              call cardin(card)
              elist=logkey(card,'no-energy-list',.false.,' ')
              nen=intkey(card,'number-of-energies',1,' ')
              if (.not.elist) then
                  call fparr(card,'energies',energy,nen,' ')
              else
                  e0=fpkey(card,'first-energy',1.e-05,' ')
                  deltae=fpkey(card,'energy-step',.1d0,' ')
                  do 510 i=1,nen
                     energy(i)=e0
                     e0=e0+deltae
 510              continue   
              endif  
              dopse=.true.
          endif
          do 600 i=1,nen
             call rmtrx(z(rmat),z(gamic),z(eig),energy(i),nc,newmat,
     1                  prntld)
             call match(z(rmat),z(kmat),z(coef),z(coef),z(tmtrx),
     1                  z(cross),z(y0),z(dy0),z(y1),z(dy1),z(jbox),
     2                  z(jpbox),z(ybox),z(ypbox),z(wronbx),z(scr1),
     3                  ia(ipvt),lval,ec,energy(i),rbox,nc,nopen,
     4                  prntkm)
             if (nopen.ne.0) then
                 call egnpse(z(kmat),z(coef),psesum,z(dum),nc,nopen,
     1                       energy(i),prntep)
             endif
  600     continue
      endif     
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)               
      stop
    1 format(/,5x,'integration region = ',i3,1x,'number points = ',i4,
     1       /,5x,'starting value     = ',e15.8,1x,'ending value = ',
     2             e15.8)
    2 format(/,5x,'   point  ',6x,'    value    ',10x,'    weight    ')
    3 format(8x,i3,8x,e15.8,8x,e15.8)
    4 format (/5x,'channel information for 'i4,
     1            ' coupled channel calculation.',/,20x,
     2            'R-matrix radius = ', e15.8)
c    5 format (/,5x,'chnl',2x,'l',2x,'zero val. bnfs.',2x,
c     1          'zero deriv. bfns.',2x,'exp. bfns.',6x,  
c     2          'energy'/,'     ----  -  ---------------',
c     3                    '  -----------------  ----------      ------')
    5 format (/,5x,'chnl',2x,'l',2x,'zero val. bnfs.',2x,
     1          'zero deriv. bfns.',6x,  
     2          'energy'/,'     ----  -  ---------------',
     3                    '  -----------------      ------')
    6 format(5x,i3,2x,i2,7x,i3,13x,i3,10x,e15.8)
c    6 format(5x,i3,2x,i2,7x,i3,13x,i3,14x,i3,6x,e15.8)
    7 format(/,'     sum of integration weights = ',e15.8)
    8 format(/,25x,'Summary of Scattering Calculation',
     1      /,'    energy    ',8x,'eigenphase sum' )
    9 format(e15.8,5x,e15.8)           
      end





