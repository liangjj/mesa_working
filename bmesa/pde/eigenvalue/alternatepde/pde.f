*deck pde.f 
c***begin prologue     pde
c***date written       970810   (yymmdd)
c***revision date               (yymmdd)
c***keywords           pde
c***                   
c***author             schneider, b. i.(nsf)
c***source             pde
c***purpose            three dimensional eigenvalue code
c
c***description        the one, two or three dimensional eigenvalue
c***                   problem is solved using a discrete variable
c***                   representation. 
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       pde
      program pde
c
      implicit integer (a-z)
      parameter ( maxgrd=20 )
      character*4096 ops
      character*8 prtflg
      character*2 itoc, ic, jc
      character*3 ans
      character*80 cpass, title, chrkey, prtit
      character*32 titphr, qdtyp
      character*1600 card
      character*128 filham, filbec
      character*2 atom
      character*24 mattyp
      character*16 coord
      character*8 qtyp
      logical dollar, logkey, prnh 
      logical prnht, prall, nfix, prnh0
      logical toau, useau, hamd, incore, cgrid
      logical itdiag, prbufh
      logical prdvd, dvdall, frmdsk
      logical bypass, scatt, typot, totrap
      real*8 z, y, scr, w
      real*8 amass, omegat, bigomg, left, right, energy 
      real*8 pi, dleft, dright
      real*8 thresh, cnverg, hbar
      real*8 massau, lenau, timau, dum, scale, sfac, lenscl, escal
      complex*16 cdum
      pointer(p1,z(1)),(p1,iz(1))
      pointer(p2,y(1)),(p2,iy(1))
      pointer(pscr,scr(1))
      pointer(p3,w(1)),(p3,iw(1))
      dimension nmax(3), npt(3), q(3), wt(3), p(3), dp(3), ddp(3)
      dimension zq(3), zwt(3), zp(3), zdp(3), zddp(3)
      dimension zeig(3), zham(3), zv(3), zvecl(3), zvecr(3)
      dimension omegat(3)
      dimension prnh(3)
      dimension prdvd(11), qtyp(3), qdtyp(3)
      dimension sfac(3), energy(100)
      dimension tmp(5), zscr(4)
      dimension idum(maxgrd), rdum(maxgrd,maxgrd)
      common/io/inp, iout      
      data pi/3.14159265358979323844d0/
c     hbar in joule-sec      
      data hbar/1.054592d-34/
c
c          mass in Kilograms length in meters
      data massau, lenau, timau / 9.109558d-31, 5.291771d-11, 
     1                            2.418884d-17 /
      call drum
      write(iout,*)
      write(iout,*) '    three dimensional pde eigenvalue code       '
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)
      call iosys ('read character "bec filename" from rwf',
     1             -1,0,0,filbec)
      call iosys ('open bec as old',0,0,0,filbec)
c
c                       set options
c
      numdim=intkey(ops,'number-of-dimensions',1,' ')   
      scatt=logkey(ops,'scattering-calculation',.false.,' ')
      typot=logkey(ops,'coulomb-potential',.false.,' ')
      itdiag=logkey(ops,'iterative-diagonalization',.false.,' ')      
      mattyp=chrkey(ops,'matrix-type','real-symmetric',' ')
      prall=logkey(ops,'print=m7020=all',.false.,' ')
      bypass=logkey(ops,'bypass=on',.false.,' ')
c----------------------------------------------------------------------c
c   1. if task = diagonalize-h0, the code will diagonalize a hamiltonian
c      consisting of the kinetic energy plus one and two body local 
c      interactions.  the interaction terms are not assumed to depend on
c      the wavefunction.
c----------------------------------------------------------------------c 
      write(iout,1)
      if(bypass) then
         call genmat(numdim,itdiag,mattyp)
         call chainx(0)               
         stop
      endif   
      if(prall) then
         prnh0=.true. 
         prnht=.true. 
         prnh(1)=.true.
         prnh(2)=.true.
         prnh(3)=.true.
      else         
         prnh0=logkey(ops,'print=m7020=h0',.false.,' ')
         prnht=logkey(ops,'print=m7020=h',.false.,' ')
         prnh(1)=logkey(ops,'print=m7020=q1-hamiltonian',.false.,' ')
         prnh(2)=logkey(ops,'print=m7020=q2-hamiltonian',.false.,' ')
         prnh(3)=logkey(ops,'print=m7020=q3-hamiltonian',.false.,' ')
      endif
      if(scatt) then
         nen=intkey(ops,'number-of-energies',1,' ')
         call fparr(ops,'energies',energy,nen,' ')
      endif            
      toau=logkey(ops,'to-atomic-units',.false.,' ')
      useau=logkey(ops,'use-atomic-units',.false.,' ')
      totrap=logkey(ops,'use-trap-units',.false.,' ')

c
c              set various program options
c      
c
c         set trap configuration, numbers of atoms, various constants
c         and other physical parameters
c       
      amass=1.d0
      atom=chrkey(ops,'atom','cs',' ')
      if( dollar('$trap',card,cpass,inp) ) then
          call atmdat(atom,amass,omegat)
          if(toau) then
              write(iout,*) '          converting to atomic units'
              omegat(1)=omegat(1)*timau
              omegat(2)=omegat(2)*timau
              omegat(3)=omegat(3)*timau
              amass=amass/massau
              scatl=scatl/lenau
          endif
          if(useau) then
             write(iout,*) '          assuming atomic units'
             amass=1.d0
             omegat(1)=1.d0
             omegat(2)=1.d0
             omegat(3)=1.d0
          endif
          if(totrap) then
             bigomg=max(omegat(1),omegat(2),omegat(3))
             omegat(1)=omegat(1)/bigomg
             omegat(2)=omegat(2)/bigomg
             omegat(3)=omegat(3)/bigomg
             lenscl=sqrt(hbar/(amass*bigomg))
             escal=hbar*bigomg
          endif
          if(numdim.eq.1) then
             write(iout,2) amass, omegat(1)
          elseif(numdim.eq.2) then
             write(iout,3) amass, (omegat(i),i=1,2)
          elseif(numdim.eq.3) then
             write(iout,4) amass, (omegat(i),i=1,3)
          endif
      endif
c
c               spatial basis set information
c
      words=1
      n3d=1
      maxd=0
      do 10 i=1,numdim
         ic=itoc(i)
         qtyp(i)=chrkey(ops,'grid-type-dimension-'//ic,'x',' ')
         ug=intkey(ops,'grid-number-dimension-'//ic,1,' ')
         call iosys('read integer "no. subgrids for '
     1              //qtyp(i)//'" from bec',1,nsubg,0,' ')
         call iosys('read integer "left fn value for '
     1              //qtyp(i)//'" from bec',1,fl,0,' ')
         call iosys('read integer "right fn value for '
     1              //qtyp(i)//'" from bec',1,fr,0,' ')
         call iosys('read integer "no. points for '
     1              //qtyp(i)//'" from bec',nsubg,idum,0,' ')
         npt(i)=idum(ug)
         maxd=max(maxd,npt(i))
         call iosys('read integer "mod. no. points '//
     1              'for '//qtyp(i)//'" from bec',nsubg,idum,0,' ')
         nmax(i)=idum(ug)
         titphr=qtyp(i)//' pointer'
         call iosys('read character "'//titphr//'" from bec',
     1               0,0,0,title)
         call iosys('rewind '//title//' on bec read-and-write',0,
     1               0,0,' ')
         call iosys('read integer '//title//' from bec '//
     1              'without rewinding',nsubg,idum,2*nsubg,' ')
         q(i)=idum(ug)
         call iosys('read integer '//title//' from bec without '//
     1              'rewinding',nsubg,idum,0,' ')
         wt(i)=idum(ug)
         call iosys('read integer '//title//' from bec without '//
     1              'rewinding',maxgrd*maxgrd,rdum,3*maxgrd*maxgrd,' ')
         p(i)=rdum(ug,ug)
         call iosys('read integer '//title//' from bec without '//
     1              'rewinding',maxgrd*maxgrd,rdum,0,' ')
         dp(i)=rdum(ug,ug)
         call iosys('read integer '//title//' from bec without '//
     1              'rewinding',maxgrd*maxgrd,rdum,0,' ')
         ddp(i)=rdum(ug,ug)
         n3d=n3d*nmax(i)
         zq(i)=words
         zwt(i)=zq(i)+npt(i)
         zp(i)=zwt(i)+npt(i)
         zdp(i)=zp(i)+npt(i)*npt(i)
         zddp(i)=zdp(i)+npt(i)*npt(i)
         zeig(i)=zddp(i)+npt(i)*npt(i)
         zham(i)=zeig(i)+npt(i)
         zv(i)=zham(i)+npt(i)*npt(i)
         zvecl(i)=zv(i)+npt(i)
         zvecr(i)=zvecl(i)+npt(i)*npt(i)
         words=zvecr(i)+npt(i)*npt(i)
 10   continue
      words=wpadti(words)
      call memory(words,p1,ngot1,'functions',0)
      do 20 i=1,numdim
         ic=itoc(i)
         titphr=qtyp(i)//' title'
         call iosys('read character "'//titphr//'" from bec',
     1               0,0,0,title)
         call iosys('read real '//title//' from bec',npt(i),
     1               z(zq(i)),q(i)-1,' ')
c         prtit='points'
c         call prntrm(prtit,z(zq(i)),nmax(i),1,nmax(i),1,iout)
         call iosys('read real '//title//' from bec',npt(i),
     1               z(zwt(i)),wt(i)-1,' ')
c         prtit='weights'
c         call prntrm(prtit,z(zwt(i)),nmax(i),1,nmax(i),1,iout)
         call iosys('read real '//title//' from bec',
     1               npt(i)*npt(i),z(zp(i)),p(i)-1,' ')
         call iosys('read real '//title//' from bec',
     1               npt(i)*npt(i),z(zdp(i)),dp(i)-1,' ')
         call iosys('read real '//title//' from bec',
     1               npt(i)*npt(i),z(zddp(i)),ddp(i)-1,' ')
c         prtit='p'
c         call prntrm(prtit,z(zp(i)),nmax(i),nmax(i),
c     1               nmax(i),nmax(i),iout)
c         prtit='dp'
c         call prntrm(prtit,z(zdp(i)),nmax(i),nmax(i),
c     1              nmax(i),nmax(i),iout)
c         prtit='ddp'
c         call prntrm(prtit,z(zddp(i)),nmax(i),nmax(i),
c     1               nmax(i),nmax(i),iout)
         call iosys('read real "ke eigenvalues-'//ic//' for '
     1              //qtyp(i)//'" from bec',nmax(i),z(zeig(i)),0,' ')
         prtit='eigenvalues'
         call prntrm(prtit,z(zeig(i)),nmax(i),1,nmax(i),1,iout)
         call iosys('read real "v mtrx-'//ic//' for '
     1              //qtyp(i)//'" from bec',nmax(i),z(zv(i)),0,' ')
         call iosys('read real "ke mtrx-'//ic//' for '//qtyp(i)//
     1              '" from bec',nmax(i)*nmax(i),z(zham(i)),0,' ')           
         call iosys('read real "trn mtrx-'//ic//' for '//qtyp(i)//
     1              '" from bec',nmax(i)*nmax(i),z(zvecl(i)),0,' ')
         call mkone(z(zham(i)),z(zv(i)),nmax(i))
         call copy(z(zvecl(i)),z(zvecr(i)),nmax(i)*nmax(i))
 20   continue
      vtot=1
      words=wpadti(vtot+n3d)
      call memory(words,p2,ngot2,'pot',0)            
      words=wpadti(1+n3d*numdim)
      call memory(words,pscr,nscr,'scr',0)            
      call vmat(z(zq(1)),z(zq(2)),z(zq(3)),y(vtot),scr(1),hbar,amass,
     1          scale,n3d,nmax(1),nmax(2),nmax(3),numdim,prnh0,ops)
      call memory(-nscr,pscr,junk,'scr',junk)   

c
c     subtract the one body interactions included in the zeroth
c     order hamiltonian.
c
      call vpert(y(vtot),z(zv(1)),z(zv(2)),z(zv(3)),n3d,
     1           nmax(1),nmax(2),nmax(3),numdim,mattyp)
      nroots=n3d         
c
c             set diagonalization procedures
c      
      if(itdiag) then
         call dvddat(card,cpass,n3d,nroots,ntrials,nattim,cnverg,
     1               thresh,niter,nvec,lenbuf,cgrid,prbufh,hamd,
     2               n0,prdvd,dvdall,filham)
         write(iout,5) nroots, nattim, thresh, cnverg, niter
      else
         nroots=intkey(ops,'number-of-roots',n3d,' ')
         nroots=min(nroots,n3d)
         write(iout,6) n3d, nroots
      endif
      zscr(1)=1
      fac=1
      lwork=0
      mwork=0
      if(mattyp.eq.'complex'.or.
     1          mattyp.eq.'real-unsymmetric') then
         fac=2
         lwork=20*maxd
         mwork=10*nvec*2
      endif
      zscr(2)=zscr(1)+max(2*fac*maxd*maxd,3*fac*n3d,lwork)
      words=zscr(2)+2*fac*maxd*maxd
      if(itdiag) then
         zscr(3)=words
         zscr(4)=zscr(3)+max(n3d*fac*nvec,mwork)
         words=zscr(4)+max(n3d*fac*nvec,2*fac*nvec)
      endif
      words=wpadti(words)
      call memory(words,pscr,nscr,'scr',0)
      words=1
      if(.not.itdiag) then
         eigtot=words
         hamtot=eigtot+fac*n3d
         bigvec=hamtot+fac*n3d*n3d
         bigvec2=bigvec+fac*n3d*n3d
         words=bigvec2+fac*n3d*n3d
      else
         hbuf=words
         ibuf=wpadti(hbuf+fac*lenbuf)
         diag=iadtwp(ibuf+2*lenbuf)
         etrial=diag+n3d*fac
         trials=etrial+n3d*fac
         diagv=trials+n3d*ntrials*fac
         psi=diagv+n3d*fac
         pvec=psi+n3d*nroots*fac
         hpvec=pvec+n3d*nvec*fac
         vec=hpvec+nvec*n3d*fac
         bmat=vec+n3d*nvec*fac
         bmatm=bmat+nvec*nvec*fac
         svec=bmatm+nvec*nvec*fac
         sveca=svec+nvec*nvec*fac
         eigtot=sveca+nvec*nvec*fac
         etmp=eigtot+n3d*fac
         words=etmp+n3d*fac
      endif      
      ind=wpadti(words)
      itmp=ind+(numdim+1)*n3d
      words=itmp+numdim*n3d            
      call memory(words,p3,ngot3,'diag',0)
      lwork=max(lwork,mwork)
      call indset(iw(ind),nmax,n3d,numdim)
      if(itdiag) then
         call iosys('write integer "length of davidson vector" '//
     1              'to ham',1,n3d,0,' ')
         call vtrial(z(zham(1)),z(zham(1)),
     1               z(zham(2)),z(zham(2)),
     2               z(zham(3)),z(zham(3)),
     3               z(zvecl(1)),z(zvecl(1)),
     4               z(zvecl(2)),z(zvecl(2)),
     5               z(zvecl(3)),z(zvecl(3)),
     6               z(zeig(1)),z(zeig(1)),
     7               z(zeig(2)),z(zeig(2)),
     8               z(zeig(3)),z(zeig(3)),
     9               w(trials),w(trials),
     x               w(etrial),w(etrial),
     x               iw(itmp),nmax(1),nmax(2),nmax(3),
     x               numdim,n3d,ntrials,.false.,mattyp)         
         title='buffered hamiltonian for '//itoc(numdim)
     1         //'d dimensional hamiltonian'
         call hamil(z(zham(1)),z(zham(1)),
     1              z(zham(2)),z(zham(2)),
     2              z(zham(3)),z(zham(3)),
     3              y(vtot),w(hbuf),w(hbuf),iw(ibuf),
     4              iw(ind),w(diag),w(diag),lenbuf,
     5              nmax(1),nmax(2),nmax(3),numdim,n3d,
     6              prbufh,nel,incore,title,mattyp)               
         if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
            call cdvd(w(hbuf),iw(ibuf),w(diag),w(trials),
     1                z(zeig(1)),z(zeig(2)),z(zeig(2)),
     2                z(zvecl(1)),z(zvecl(2)),z(zvecl(3)),
     3                z(zvecr(1)),z(zvecr(2)),z(zvecr(3)),
     4                w(pvec),w(hpvec),w(vec),
     5                w(bmat),w(bmatm),w(etmp),w(svec),w(sveca),
     6                w(vec),w(eigtot),scr(zscr(3)),scr(zscr(4)),
     7                lwork,scale,cnverg,thresh,
     8                nmax(1),nmax(2),nmax(3),
     9                numdim,n3d,nroots,ntrials,niter,nvec,lenbuf,
     x                hamd,incore,nel,prdvd,mattyp)
         elseif(mattyp.eq.'real-symmetric') then
            call rsdvd(w(hbuf),iw(ibuf),w(diag),w(trials),
     1                 z(zeig(1)),z(zeig(2)),z(zeig(2)),
     2                 z(zvecl(1)),z(zvecl(2)),z(zvecl(3)),w(psi),
     3                 w(pvec),w(hpvec),w(vec),
     4                 w(bmat),w(bmatm),w(etmp),w(svec),w(vec),
     5                 w(eigtot),scr(zscr(3)),scr(zscr(4)),
     6                 scale,cnverg,thresh,
     7                 nmax(1),nmax(2),nmax(3),
     7                 numdim,n3d,nroots,ntrials,nattim,niter,
     8                 nvec,lenbuf,hamd,incore,nel,prdvd,mattyp)
         endif          
      else    
         call lschr(z(zham(1)),z(zham(1)),
     1              z(zham(2)),z(zham(2)),
     2              z(zham(3)),z(zham(3)),
     3              y(vtot),w(hamtot),w(hamtot),
     4              w(eigtot),w(eigtot),
     5              w(bigvec),w(bigvec),w(bigvec2),
     6              w(bigvec),w(bigvec2),
     7              scr(zscr(1)),lwork,scr(zscr(2)),
     8              iw(ind),n3d,nmax(1),nmax(2),nmax(3),
     9              numdim,.true.,.true.,prnht,mattyp)                      
      endif
      if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
         call cvscal(w(eigtot),w(eigtot),scale,nroots)
         title='eigenvalues of hamiltonian'
         call prntcm(title,w(eigtot),nroots,1,nroots,1,iout)
      elseif(mattyp.eq.'real-symmetric') then
         call vscale(w(eigtot),w(eigtot),scale,nroots)
         title='eigenvalues of hamiltonian'
         call prntfm(title,w(eigtot),nroots,1,nroots,1,iout)
      endif
      if(scatt) then
         do 200 i=1,nen
            call rmtrx(z(zp(1)),z(zq(1)),w(eigtot),
     1                 w(bigvec),energy(i),typot,n3d)
  200    continue         
      endif
      call memory(-ngot1,p1,junk,'functions',junk)
      call memory(-ngot2,p2,junk,'pot',junk)            
      call memory(-nscr,pscr,junk,'scr',junk)
      call memory(-ngot3,p3,junk,'diag',junk)
      call chainx(0)               
      stop
 1    format(/,15x,'program options',/,/,5x,
     1             'diagonalize zeroth-order hamiltonian')  
 2    format(/,15x,'trap configuration and atomic parameters',/,/,5x,
     1             'atomic mass                      = ',e15.8,/,5x,
     2             'trap q1 frequency                = ',e15.8)
 3    format(/,15x,'trap configuration and atomic parameters',/,/,5x,
     1             'atomic mass                       = ',e15.8,/,5x,
     2             'trap q1 frequency                 = ',e15.8,/,5x,
     3             'trap q2 frequency                 = ',e15.8)
 4    format(/,15x,'trap configuration and atomic parameters',/,/,5x,
     1             'atomic mass                       = ',e15.8,/,5x,
     2             'trap q1 frequency                 = ',e15.8,/,5x,
     3             'trap q2 frequency                 = ',e15.8,/,5x,
     4             'trap q3 frequency                 = ',e15.8)
 5    format(/,15x,'iterative diagonalization information',/,/,5x,
     1             'number of roots                    = ',i3,/,5x,
     2             'number of roots at a time          = ',i3,/,5x
     3             'overlap tolerance                  = ',e15.8,/,5x,
     4             'convergence criterion              = ',e15.8,/,5x,
     5             'maximum number of iterations       = ',i6)
 6    format(/,15x,'direct diagonalization',/,/,5x,
     1             'size of matrix  = ',i3,/,5x,
     2             'number of roots = ',i3)
      end
