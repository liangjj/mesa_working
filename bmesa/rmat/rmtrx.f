*deck rmtrx 
c***begin prologue     m6250
c***date written       960128   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6250, link 6250
c***author             schneider, b. i.(nsf)
c***source             m6250
c***purpose            coupled channel r-matrix code using orthogonal 
c***                   functions as translational basis.
c***description        orthogonal polynomials are used as a basis for the
c***                   translational scattering functions.
c***                   
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m6250
      program rmtrx
c
      implicit integer (a-z)
      parameter ( maxl=10, ptmax=1000, chnmax=20 )
      parameter ( noene=200 )
      dimension nobf(chnmax), energy(noene)
      dimension lval(chnmax), ec(chnmax), wt(ptmax), point(ptmax)
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
      logical elist, preigs, scat, dvr
      common z(1)
      dimension ia(1)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
      real*8 z, rbox, fpkey, root, droot, wt, point, rst, rfn, rfs
      real*8 ec, wtsum, e0, deltae, energy, psesum, vcij, range
      real*8 rfinal, left, right
c
      call drum
      write(iout,*)
      write(iout,*) '                         m6250:coupled channel'//
     1                                        ' r-matrix link'
      call iosys ('read character options from rwf',-1,0,0,ops)
      prntwv=logkey(ops,'print=m6250=wavefunction',.false.,' ')
      prntgr=logkey(ops,'print=m6250=grid',.false.,' ')
      prntov=logkey(ops,'print=m6250=overlap-matrix',.false.,' ')
      prntke=logkey(ops,'print=m6250=kinetic-energy-matrix',.false.,' ')
      prntpe=logkey(ops,'print=m6250=potential-energy-matrix',
     1              .false.,' ')
      preigs=logkey(ops,'print=m6250=hamiltonian-eigenvalues',
     1                   .false.,' ')
      prnth=logkey(ops,'print=m6250=hamiltonian-matrix',.false.,' ')
      prntld=logkey(ops,'print=m6250=r-matrix',.false.,' ')
      prntkm=logkey(ops,'print=m6250=k-matrix',.false.,' ')
      prntep=logkey(ops,'print=m6250=eigenphase-vectors',.false.,' ')
      prntir=logkey(ops,'print=m6250=intermediate-results',
     1              .false.,' ')
      pottyp=chrkey(ops,'potential-type','square-well',' ')
      scat=logkey(ops,'m6250=scattering=on',.false.,' ')
      dvr=logkey(ops,'m6250=discrete-variable-representation',.
     1           false.,' ')
      call iosys ('read character "linear algebraic filename" from rwf',
     1            -1,0,0,fillam)
      call iosys ('open lamdat as old',0,0,0,fillam)
      if (posinp('$grid',cpass) ) then
          call cardin(card)
          rbox=fpkey(card,'r-matrix-radius',10.d0,' ')
          rfinal=fpkey(card,'final-r',rbox,' ')
          range=fpkey(card,'potential-range',rbox,' ')
      endif
      call iosys ('write character "potential type" to lamdat',0,0,
     1             0,pottyp)
      call iosys('write real "r-matrix radius" to lamdat',1,rbox,0,' ')
      call iosys('write real "final radius" to lamdat',1,rfinal,0,' ')
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
      ij=0
      mxchn=0
      matsiz=0
      do 30 i=1,nc
         mxchn=max(mxchn,nobf(i))
         matsiz=matsiz+nobf(i)
         call pakstr(itoc(i),ii)
         do 40 j=1,i
            ij=ij+1
            call pakstr(itoc(j),jj)
            chr1=itoc(i)
            chr2=itoc(j)
            vcij(ij)=fpkey(card,'v'//chr1(1:ii)//chr2(1:jj),-1.d0,' ')
 40      continue
 30   continue
      write(iout,4) nc, rbox
      write(iout,5)
      lmax=0
      do 50 i=1,nc
         lmax=max(lmax,lval(i))
         ii=i
         write(iout,6) i, lval(i), nobf(i), ec(i)
 50   continue
      call iosys ('write integer "maximum l value" to lamdat',1,
     1             lmax,0,' ') 
      call iosys('read integer "number of polynomials" from lamdat',
     1            1,npoly,0,' ')
      call iosys('read real "left boundary" from lamdat',1,left,0,' ')
      call iosys('read real "right boundary" from lamdat',1,right,0,' ')    
      lp=lmax+1     
      ioff=1
      do 60 i=1,2
         f=ioff
         df=f+ptcnt*matsiz
         ddf=df+ptcnt*matsiz
         a=ddf+ptcnt*matsiz
         b=a+npoly+1
         flst=b+npoly+1
         dflst=flst+matsiz 
         x=dflst+matsiz 
         scr=x+ptcnt
         w1=scr+ptcnt
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
   60 continue
c      calculate the channel basis functions
      f0=f
      df0=df
      ddf0=ddf
      do 70 i=1,nc
         call mkcfun(z(f0),z(df0),z(ddf0),point,z(a),z(b),left,right,
     1               lval(i),nobf(i),ptcnt,npoly,prntwv)
         f0=f0+ptcnt*nobf(i)
         df0=df0+ptcnt*nobf(i)
         ddf0=ddf0+ptcnt*nobf(i)
 70   continue
c      call lnkerr('quit')   
      call prep(z(f),z(df),z(ddf),z(flst),z(dflst),z(scr),wt,ptcnt,
     1          matsiz)
c      calculate the overlap and kinetic energy matrix elements.
c      they are diagonal in the channel index.
      f0=f
      df0=df
      ddf0=ddf
      flst0=flst
      dflst0=dflst
      t0=t
      do 300 i=1,nc
c        calculate the overlap matrix of the primitive functions.
c        in principle this is already a diagonal set.
         call smat(z(f0),z(s),ptcnt,nobf(i),prntir)
c        normalize the overlap
         call renrm(z(s),z(f0),z(df0),z(ddf0),z(flst0),z(dflst0),
     1              z(scr1),ptcnt,nobf(i),prntir)
c        calculate the kinetic energy matrix     
         call tmat(z(f0),z(df0),z(ddf0),z(flst0),z(dflst0),z(t0),ptcnt,
     1             nobf(i),prntir)
         if (prntov) then
             title='overlap matrix for channel = '//itoc(i)
             call prntrm(title,z(s),nobf(i),nobf(i),nobf(i),
     1                   nobf(i),iout)
         endif
         if (prntke) then
             title='kinetic energy matrix for channel = '//itoc(i)
             call prntrm(title,z(t0),nobf(i),nobf(i),nobf(i),
     1                   nobf(i),iout)
         endif
         f0=f0+ptcnt*nobf(i)
         df0=df0+ptcnt*nobf(i)
         ddf0=ddf0+ptcnt*nobf(i)
         flst0=flst0+nobf(i)
         dflst0=dflst0+nobf(i)
         t0=t0+nobf(i)*nobf(i)
  300 continue
      fi=f
      vij=v
      ntri=nc*(nc+1)/2
      do 400 i=1,nc
         fj=f
         do 500 j=1,i
            call potntl(z(fi),z(fj),vcij,z(pot),z(vij),point,range,
     1                  nobf(i),nobf(j),ptcnt,pottyp,i,j,ntri,prntir)
            if (prntpe) then
                call pakstr(itoc(i),ii)
                call pakstr(itoc(j),jj)
                chr1=itoc(i)
                chr2=itoc(j)
                title='potential energy matrix for channels = ('
     1                //chr1(1:ii)//','//chr2(1:jj)//')'
                call prntrm(title,z(vij),nobf(i),nobf(j),nobf(i),
     1                      nobf(j),iout)
            endif                       
            fj=fj+nobf(j)*ptcnt
            vij=vij+nobf(i)*nobf(j)
  500    continue
         fi=fi+nobf(i)*ptcnt
  400 continue
      call filham(z(ham),z(t),z(v),ec,nobf,nc,matsiz)
      if (prnth) then
          title='final hamiltonian matrix'
          call prntrm(title,z(ham),matsiz,matsiz,matsiz,matsiz,iout)
      endif           
      call diag(z(ham),z(eig),z(dum),matsiz,preigs)
      call surpsi(z(ham),z(flst),z(gamic),nobf,nc,matsiz,prntir)
      dopse=.false.
      if ( scat ) then
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
             call conrmt(z(rmat),z(gamic),z(eig),energy(i),nc,newmat,
     1                   prntld)
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
    5 format (/,5x,'chnl',2x,'l',2x,'no. of bsis fns',6x,  
     1          'energy'/,'     ----  -  --------------',
     2                    '      ------')
    6 format(5x,i3,2x,i2,7x,i3,10x,e15.8)
    7 format(/,'     sum of integration weights = ',e15.8)
    8 format(/,25x,'Summary of Scattering Calculation',
     1      /,'    energy    ',8x,'eigenphase sum' )
    9 format(e15.8,5x,e15.8)           
      end





