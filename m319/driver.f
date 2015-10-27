*deck @(#)driver.f	5.3  11/28/95
      subroutine driver(ptprim,noprim,nbtype,nocont,ptcont,ex,cont,
     $                  ncont,c,ints,labels,lenbuf,start,nx,ny,nz,
     $                  lenxyz,nobf,nocart,mintyp,minmom,maxmom,
     $                  z,maxcor,nat,nbasis,nnp,nprim,ops,cutexp,prnt,
     $                  nintgr,dis,atnam,drop,pkindx)
c
      implicit integer(a-z)
c
      character*(*) ops, atnam(nat)
      real*8 ex(nprim),cont(ncont),c(3,nat),ints(lenbuf),z(maxcor)
      real*8 cutexp,pi252,pcnt1,pcnt2
      real*8 one,uh,uhe,up,v,dt,exchng,fpkey,rmax,dis(nat,nat)
      real*8 ud,deltat
      real*8 zero
      integer ptprim(nat,nbtype),noprim(nat,nbtype),nocont(nat,nbtype)
      integer ptcont(nat,nbtype),labels(nintgr,lenbuf),nobf(nbtype)
      integer start(nat,nbtype),nocart(nbtype),mintyp(nbtype)
      integer minmom(nbtype),maxmom(nbtype)
      integer nx(lenxyz),ny(lenxyz),nz(lenxyz)
      integer symcen(4),angmom(4)
      integer pkindx(nbasis)
      logical prnt,debug,logkey
      logical drop
c
      common /io/     inp,iout
      parameter (debug=.false.)
c
c     ----- get 2*pi**(5/2) -----
c
      call iosys('read real pi from rwf',1,pi252,0,' ')
      pi252=2.0d+00*pi252**2.5d+00
      zero=0.0d+00
      one=1.0d+00
c
c     get the parameters for the ppp model.
      uh=fpkey(ops,'ppp=uh',0.0d+00,' ')
      if(uh.eq.0.0d+00) uh=fpkey(ops,'ppp=u',0.0d+00,' ')
      uhe=fpkey(ops,'ppp=uhe',0.0d+00,' ')
      up=fpkey(ops,'ppp=up',0.0d+00,' ')
      v=fpkey(ops,'ppp=v',0.0d+00,' ')
      dt=fpkey(ops,'ppp=dt',0.0d+00,' ')
      exchng=fpkey(ops,'ppp=k',0.0d+00,' ')
      rmax=fpkey(ops,'ppp=rmax',0.0d+00,' ')
c
c     --- get parameters for two-site model.  this model has only 'on-site'
c         two electron integrals.  it is run by using the he atom
c         with two s-type basis functions.
      rmax=fpkey(ops,'twosit=rmax',rmax,' ')
      ud=fpkey(ops,'twosit=ud',0.0d+00,' ')
      deltat=fpkey(ops,'twosit=deltat',0.0d+00,' ')
c
c     check for inadmissable input.
      do 50 i=1,nat
         do 40 itype=1,nbtype
            nprimi=noprim(i,itype)
            if (nprimi.gt.0.and.itype.gt.1) then
               if(dt.ne.zero.or.exchng.ne.zero)
     $            call lnkerr('  (ii|ij)/(ij|ij) terms only allowed '
     $                       //'with s-functions')
            end if
   40    continue
   50 continue
c
c     ----- initialize counters -----
c
      maxc=0
      maxc1=0
      maxc2=0
      maxc3=0
      maxc4=0
      numbuf=1
      ntotal=0
      nonzer=0
      nactul=0
c
c     ----- check amount of core needed -----
c
      call       sizer1(ptprim,noprim,nbtype,nocont,ptcont,ex,cont,
     #                  ncont,c,ints,labels,lenbuf,start,nx,ny,nz,
     #                  lenxyz,nobf,nocart,mintyp,minmom,maxmom,
     #                  z,maxcor,nat,nbasis,nnp,nprim,ops,cutexp,need1)
c
      call       sizer2(ptprim,noprim,nbtype,nocont,ptcont,ex,cont,
     #                  ncont,c,ints,labels,lenbuf,start,nx,ny,nz,
     #                  lenxyz,nobf,nocart,mintyp,minmom,maxmom,
     #                  z,maxcor,nat,nbasis,nnp,nprim,ops,cutexp,need2)
c
c     ----- check that can get enough core -----
c
      call getscm(0,z,canget,'m319:driver',0)
      if(prnt) then
         minlen=need1
         maxlen=need2
         write(iout,191) minlen,maxlen
  191    format(5x,'minimum necessary field length ',i9,
     $         /5x,'maximum usable field length    ',i9)
      endif
      if (wpadti(need1).gt.canget) 
     $           call lnkerr('not enough core available for m319')
c
c     add a little just to be on the safe side.
      needed=min(canget,wpadti(need2+1000))
c
c     ----- get the memory -----
c
      call getscm(needed,z,junk,'m319:driver',0)
      maxcor=iadtwp(junk)
c
c     compute the distance matrix.
      call dismat(nat,c,dis,one)
c
      if (debug) then
         write(iout,*) ' two-electron integrals'
      end if
      do 8000 iatom=1,nat
         symcen(1)=iatom
         do 7000 itype=1,nbtype
            nprimi=noprim(iatom,itype)
            if (nprimi.le.0) go to 7000
            angmom(1)=itype
            imax=maxmom(itype)
            nconti=nocont(iatom,itype)
            do 6000 jatom=1,iatom
               symcen(2)=jatom
               if (jatom.eq.iatom) then
                  jtmax=itype
               else
                  jtmax=nbtype
               end if
               do 5000 jtype=1,jtmax
                  nprimj=noprim(jatom,jtype)
                  if (nprimj.eq.0) go to 5000
                  angmom(2)=jtype
                  jmax=maxmom(jtype)
                  ncontj=nocont(jatom,jtype)
                  do 4000 katom=1,iatom
                     symcen(3)=katom
                     if (katom.eq.iatom) then
                        ktmax=itype
                     else
                        ktmax=nbtype
                     end if
                     do 3000 ktype=1,ktmax
                        nprimk=noprim(katom,ktype)
                        if (nprimk.le.0) go to 3000
                        angmom(3)=ktype
                        kmax=maxmom(ktype)
                        ncontk=nocont(katom,ktype)
                        if (katom.eq.iatom.and.ktype.eq.itype) then
                           latmax=jatom
                        else
                           latmax=katom
                        end if
                        do 2000 latom=1,latmax
                           symcen(4)=latom
                           if (katom.eq.iatom.and.ktype.eq.itype.and.
     $                          latom.eq.jatom) then
                              ltmax=jtype
                           else if (latom.eq.katom) then
                              ltmax=ktype
                           else
                              ltmax=nbtype
                           end if
                           do 1000 ltype=1,ltmax
                              npriml=noprim(latom,ltype)
                              if (npriml.le.0) go to 1000
                              angmom(4)=ltype
                              lmax=maxmom(ltype)
                              ncontl=nocont(latom,ltype)
c
               nmax=imax+jmax
               mmax=kmax+lmax
               nroots=(nmax+mmax)/2+1
               npint=nprimi*nprimj*nprimk*npriml
               ncint=nconti*ncontj*ncontk*ncontl
               nij=nprimi*nprimj
               nkl=nprimk*npriml
               lenblk=nocart(angmom(1))*nocart(angmom(2))*
     #                nocart(angmom(3))*nocart(angmom(4))
               numint=ncint*lenblk
c
c     ----- core allocation. nb there is a great deal of overlapping --
c
               a=1
               ar=a+nij
               xyza=ar+nij
               xyzam1=xyza+nij*3
               xyzam3=xyzam1+nij*3
               b=xyzam3+nij*3
               br=b+nkl
               xyzb=br+nkl
               xyzbm1=xyzb+nkl*3
               xyzbm3=xyzbm1+nkl*3
               half=xyzbm3+nkl*3
               test=half+nconti*ncontj*nkl*lenblk
c
c              ----- work out how many kl primitive sets to calculate
c                    simultaneously
c
               if (nmax+mmax.eq.0) then
                  fac1=10+5*nroots
                  fac2=3+3*nroots
               else
                  fac1=10*nroots+1+max(3*(nmax+1)*(mmax+1)*nroots,
     #                                 9+5*nroots)
                  fac2=nroots*3*(jmax+1)*(mmax+1)*((imax+1)*(lmax+1)+
     #                                             (nmax+1))+1
               end if
               fac3=numint+ncontk*npriml
               fac4=numint*2+lenblk*4
c
               space=maxcor-test-1000
               numkl=min(nkl,(space/max(fac1,fac2))/nij)
               if (numkl.le.0) call lnkerr('internal error in driver'//
     #                       ' with numkl')
               npass=(nkl+numkl-1)/numkl
               nijkl=nij*numkl
               lenv=nijkl*nroots
c
               g=test+nijkl
               ab=g
               aplusb=ab+nijkl
               urho=aplusb+nijkl
               wt=urho+lenv
               denom=wt+lenv
               t1=denom+lenv
               t2=t1+nijkl
               t3=t2+nijkl
               t4=t3+nijkl
               t5=t4+lenv
               t6=t5+lenv
               t7=t6+nijkl
               t8=t7+nijkl
               t9=t8+nijkl
c
               if (nmax+mmax.eq.0) then
                  f00=denom
                  b00=1
                  b10=1
                  bp01=1
                  c00=1
                  cp00=1
                  top1=t9+nijkl
               else
                  f00=max(g+3*(nmax+1)*(mmax+1)*lenv,t9+nijkl)
                  b00=f00+lenv
                  b10=b00+lenv
                  bp01=b10+lenv
                  c00=bp01+lenv
                  cp00=c00+lenv*3
                  top1=cp00+lenv*3
               end if
c
               if (top1.ne.test+fac1*nijkl) call lnkerr('internal '//
     #                   'error in driver with top1')
c
               if (nmax+mmax.eq.0) then
                  temp1=ab
                  i2=1
                  h=1
                  temp=1
                  tp1=1
                  top2=f00+lenv
               else
                  i2=g
                  h=i2+3*(imax+1)*(jmax+1)*(mmax+1)*(lmax+1)*lenv
                  lentmp=min(nijkl*nocart(angmom(3))*nocart(angmom(4)),
     #                        (maxcor-h+1)/2)
                  if (lentmp.lt.nijkl) call lnkerr('internal error '//
     #                 'with length of temps for vfmint')
                  temp=h
                  temp1=temp+lentmp
c                 temp=h
c                 temp1=temp+nijkl
c                 tp1=temp1+nijkl
                  top2=max(h+3*(nmax+1)*(jmax+1)*(mmax+1)*lenv,
     #                    temp1+lentmp)
               end if
c
c              if (top2.ne.test+fac2*nijkl) call lnkerr('internal '//
c    #                   ' error in driver with top2')
c
               conint=half+nconti*ncontj*nkl*lenblk
               tmp1=conint+ncint*lenblk
               lent=maxcor-tmp1+1
               if (lent.lt.ncontk*npriml) call lnkerr('internal '//
     #               'error with temporary array for trans')
               top3=tmp1+ncontk*npriml
               if (top3.ne.test+fac3) call lnkerr('internal error'//
     #                   ' in driver with top3')
c
               numint=ncint*lenblk
               plbls=half
               flbls=plbls+numint*nintgr
c               if (flbls.gt.conint) call lnkerr('internal error in '//
c     #                                  'driver with flbls')
               if (flbls+lenblk*4.le.conint) then
                  top4=conint+numint
               else
                  flbls=conint+numint
                  top4=flbls+lenblk*4
               end if
c              if (top4.ne.test+fac4) call lnkerr('internal error'//
c    #                   ' in driver with top4')
c
               top=max(top1,top2,top3,top4)
               maxc=max(top,maxc)
               maxc1=max(maxc1,top1)
               maxc2=max(maxc2,top2)
               maxc3=max(maxc3,top3)
               maxc4=max(maxc4,top4)
               if (top.gt.maxcor) call lnkerr('driver core')
c
c              ----- zero half-transformed integrals -----
c
               call rzero(z(half),nconti*ncontj*nkl*lenblk)
c
c
c              ----- loop over sets of kl primitives -----
c
               maxkl=0
               do 1001 pass=1,npass
                  minkl=maxkl+1
                  if (minkl.gt.nkl) go to 1002
                  maxkl=min(nkl,maxkl+numkl)
c
                  numkl=maxkl-minkl+1
                  nijkl=nij*numkl
                  lenv=nijkl*nroots
c
c
c
c
                   ntotal=ntotal+nijkl*lenblk
                   nbatch=0
                   nonzer=nonzer+nbatch*lenblk
c
         if(logkey(ops,'ppp',.false.,' ')) then
c           this is a fix for semiempirical ppp-like hamiltonians.
            call rzero(z(conint),numint)
            if(iatom.eq.jatom) then
               if(katom.eq.latom) then
                  if(katom.eq.iatom) then
c                    single center terms.
c                    (ii|ii)
                     if(atnam(iatom)(1:2).eq.'he') then
                        ijkl=0
                        do 25 i=1,nocart(angmom(1))
                           do 24 j=1,nocart(angmom(2))
                              do 23 k=1,nocart(angmom(3))
                                 do 22 l=1,nocart(angmom(4))
                                    ijkl=ijkl+1
                                    if((i.eq.j).and.(j.eq.k).
     $                                  and.(k.eq.l))
     $                                  z(conint+ijkl-1)=uhe
   22                            continue
   23                         continue
   24                      continue
   25                   continue
                     else if(atnam(iatom)(1:1).eq.'h') then
                        z(conint)=uh
c                     else if(atnam(iatom)(1:1).eq.'c') then
c                        ijkl=0
c                        do 25 i=1,nocart(angmom(1))
c                           do 24 j=1,nocart(angmom(2))
c                              do 23 k=1,nocart(angmom(3))
c                                 do 22 l=1,nocart(angmom(4))
c                                    ijkl=ijkl+1
c                                    if((i.eq.j).and.(j.eq.k).
c     $                                  and.(k.eq.l))
c     $                                  z(conint+ijkl-1)=up
c   22                            continue
c   23                         continue
c   24                      continue
c   25                   continue
                     else
                        call lnkerr('ppp only recognizes atoms h,he')
                     end if
                  else
c                    two-center terms.
c                    (ii|kk)
                     if(dis(iatom,katom).le.rmax) then
                        ijkl=0
                        do 35 i=1,nocart(angmom(1))
                           do 34 j=1,nocart(angmom(2))
                              do 33 k=1,nocart(angmom(3))
                                 do 32 l=1,nocart(angmom(4))
                                    ijkl=ijkl+1
                                    if((i.eq.j).and.(k.eq.l))
     $                                 z(conint+ijkl-1)=v
   32                            continue
   33                         continue
   34                      continue
   35                   continue
                     end if
                  end if
               else
                  if(katom.eq.iatom) then
c                    (ii|ij)
                     if(dis(katom,latom).le.rmax) z(conint)=dt
                  end if
               end if
            else
               if(dis(iatom,jatom).le.rmax) then
                  if(katom.eq.iatom) then
c                    (ij|ii) and (ij|ij)
                     if(latom.eq.iatom) z(conint)=dt
                     if(latom.eq.jatom) z(conint)=exchng
                  else if(katom.eq.jatom) then
c                    (ij|jj)
                     if(katom.eq.latom) z(conint)=dt
                  end if
               end if
            end if
         endif
c
c        --- this is the construction for the two-site model
c            it is invoked by the keyword twosit and specifying
c            the atoms as the lattice with two contracted 
c            s-functions per site
         if(logkey(ops,'twosit',.false.,' ')) then
            call rzero(z(conint),numint)
            if((iatom.eq.jatom).and.(katom.eq.latom)) then
               if(katom.eq.iatom) then
c                 one-center terms. first orbital is d-like
c                 (ii|ii) -- orbitals (11|11)
                  z(conint)=ud
c                 --- orbitals (21|11),(12|11),(11|21), and (11|12)
c                     hole notation  (dd|dp)
                  z(conint+1)=deltat
                  z(conint+2)=deltat
                  z(conint+4)=deltat
                  z(conint+8)=deltat
c                 electron notation (pp|pd) 
c                                   (22|21),(22|12),(21|22),(12|22)
c                 z(conint+7)=deltat
c                 z(conint+11)=deltat
c                 z(conint+13)=deltat
c                 z(conint+14)=deltat
               endif
            endif
         endif
c
               if(debug) then
                  write(iout,*) iatom,jatom,katom,latom
                  write(iout,*) (z(conint+i-1),i=1,lenblk)
               end if
c
 1001          continue
 1002          continue
c
c     ----- write the contracted integrals to tape44 -----
c
               if (nintgr.eq.1) then
c
c                 ----- 64 bit machines -----
c
                  call out64(z(conint),nconti,ncontj,ncontk,ncontl,
     #                       lenblk,symcen,angmom,start,nat,nbtype,nnp,
     #                       ints,labels,lenbuf,numbuf,
     #                       ncint,z(flbls),numint,z(plbls),
     #                       nbasis,nobf,nactul,drop,pkindx)
               else
c
c                 ----- 32 bit machines -----
c
                  call out32(z(conint),nconti,ncontj,ncontk,ncontl,
     #                       lenblk,symcen,angmom,start,nat,nbtype,nnp,
     #                       ints,labels,lenbuf,numbuf,
     #                       ncint,z(flbls),numint,z(plbls),
     #                       nbasis,nobf,nactul,drop,pkindx)
               end if
 1000                      continue
 2000                   continue
 3000                continue
 4000             continue
 5000          continue
 6000       continue
 7000    continue
 8000 continue
c
c    ----- flush final integral buffer -----
c
      nactul=nactul+numbuf-1
      ntotc=nnp*(nnp+1)/2
      ints(1)=-numbuf
      call iosys('write real "unsorted ao integrals" on rints '//
     $     'without rewinding',lenbuf,labels,0,' ')
      call iosys('write real "unsorted ao integrals" on rints '//
     $     'without rewinding',lenbuf,ints,0,' ')
c
      call iosys('endfile "unsorted ao integrals" on rints',0,0,0,' ')
cpws
      if(prnt) then
         pcnt1=float(nonzer)*100.d+00/float(ntotal)
         pcnt2=float(nactul)*100.d+00/float(ntotc)
c        write(iout,956) maxc,maxc1,maxc2,maxc3,maxc4
  956    format(5x,'maxc=',5i10)
         write(iout,957) ntotal,nonzer,pcnt1,ntotc,nactul,pcnt2
  957    format(5x,'#primitive integrals possible',1x,i10,
     $         /5x,'#primitive integrals computed',1x,i10,'(',f5.1,'%)',
     $         /5x,'#contracted integrals possible',i10,
     $         /5x,'#contracted integrals kept',4x,i10,'(',f5.1,'%)')
      endif
c
c     ----- timing -----
c
c
c
      return
      end
