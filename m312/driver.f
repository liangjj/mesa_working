*deck @(#)driver.f	5.3   4/17/95
      subroutine driver(ptprim,noprim,nbtype,nocont,ptcont,ex,cont,
     $                  ncont,c,ints,labels,lenbuf,start,nx,ny,nz,
     $                  lenxyz,nobf,nocart,mintyp,minmom,maxmom,
     $                  nat,nbasis,nnp,nprim,ops,cutexp,prnt,
     $                  nintgr,drop,pkindx)
c***begin prologue     driver.f
c***date written       yymmdd  
c***revision date      4/17/95      
c   2 july 1991        rlm at lanl
c      adding drop option
c
c***keywords           
c***author             saxe, paul(lanl)
c***source             @(#)driver.f	5.3   4/17/95
c***purpose            driver for two-electron integrals 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       driver.f
      implicit none
c     --- input variables -----
      integer maxcor,nat,nbtype,nintgr,lenxyz,nprim,ncont,ngot,idum
      integer lenbuf,nbasis,nnp
      logical drop,prnt
      character*(*) ops
      real*8 cutexp
c     --- input arrays (unmodified) ---
      integer ptprim(nat,nbtype),noprim(nat,nbtype),nocont(nat,nbtype)
      integer ptcont(nat,nbtype),labels(nintgr,lenbuf),nobf(nbtype)
      integer start(nat,nbtype),nocart(nbtype),mintyp(nbtype)
      integer minmom(nbtype),maxmom(nbtype)
      integer nx(lenxyz),ny(lenxyz),nz(lenxyz)
      integer pkindx(nbasis)
      real*8 ex(nprim),cont(ncont),c(3,nat),ints(lenbuf)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
      real*8 z
c     --- local variables ---
      integer inp,iout
      integer symcen(4),angmom(4)
      integer maxc,maxc1,maxc2,maxc3,maxc4
      integer numbuf,ntotal,nonzer,nactul
      integer need1,need2,canget,iadtwp,wpadti,wptoin
      integer minlen,maxlen,needed
      integer iatom,itype,nprimi,imax,nconti
      integer jatom,jtmax,jtype,nprimj,jmax,ncontj
      integer katom,ktmax,nprimk,ncontk,ktype,kmax
      integer latom,ltmax,ltype,npriml,ncontl,latmax,lmax
      integer nmax,mmax,lenblk,nroots,fac1,fac2,fac3,fac4
      integer a,xyza,ar,nij,br,b,nkl
      integer space,test,numkl
      integer ab,h,g,urho,aplusb,nijkl,wt,lenv,f00,b00,b10,bp01,c00,cp00
      integer i2,junk,temp,temp1,tp1,lentmp,top1,top2,top3
      integer lent,plbls,top4,conint,numint
      integer xyzam1,xyzam3,xyzb,xyzbm1,xyzbm3,maxkl,npass,minkl
      integer t1,t2,t3,t4,t5,t6,t7,t8,t9
      integer nbatch,rbatch,npint,ncint,flbls,ntotc
      integer denom,pass,half,tmp1,top
      real*8 pi252,pcnt1,pcnt2
      pointer (p,z(1))
c
      common /io/     inp,iout
c
 9000 format(5x,'minimum necessary field length ',i9,
     $      /5x,'maximum usable field length    ',i9)
 9010 format(5x,'maxc=',5i10)
 9020 format(5x,'#primitive integrals possible',1x,i10,
     $      /5x,'#primitive integrals computed',1x,i10,'(',f5.1,'%)',
     $      /5x,'#contracted integrals possible',i10,
     $      /5x,'#contracted integrals kept',4x,i10,'(',f5.1,'%)')
 9030 format(5x, 'integer words gotten = ',i10)
c
      call iosys('read integer mxcore from rwf',1,maxcor,0,' ')
      canget=iadtwp(maxcor)
c
c     --- get 2*pi**(5/2) ---
      call iosys('read real pi from rwf',1,pi252,0,' ')
      pi252=2.0d+00*pi252**2.5d+00
c
c     --- initialize counters ---
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
c     --- check amount of core needed ---
      call       sizer1(ptprim,noprim,nbtype,nocont,ptcont,ex,cont,
     $                  ncont,c,ints,labels,lenbuf,start,nx,ny,nz,
     $                  lenxyz,nobf,nocart,mintyp,minmom,maxmom,
     $                  z,canget,nat,nbasis,nnp,nprim,
     $                  ops,cutexp,need1)

      call       sizer2(ptprim,noprim,nbtype,nocont,ptcont,ex,cont,
     $                  ncont,c,ints,labels,lenbuf,start,nx,ny,nz,
     $                  lenxyz,nobf,nocart,mintyp,minmom,maxmom,
     $                  z,canget,nat,nbasis,nnp,nprim,o
     $                  ps,cutexp,need2)

c     --- check that can get enough core ---
      if(prnt) then
         minlen=need1
         maxlen=need2
         write(iout,9000) minlen,maxlen
      endif
      if (wpadti(need1).gt.maxcor)  then
         write(iout,*)'need:',need1,' real*8 words',
     $                ' can get:',canget
         call lnkerr('not enough core available'//
     $        'for m312')
      endif
c
      needed=min(canget,need2)
      needed=wptoin(needed) 
c
c     --- get the memory ---
      call getmem(needed,p,ngot,'driver',0)
      write(iout,9030) ngot
      maxcor=iadtwp(ngot)
c
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
     $                nocart(angmom(3))*nocart(angmom(4))
               numint=ncint*lenblk
c
c              --- core allocation. 
c                  nb there is a great deal of overlapping --
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
c              --- work out how many kl primitive sets to calculate
c                  simultaneously
               if (nmax+mmax.eq.0) then
                  fac1=10+5*nroots
                  fac2=3+3*nroots
               else
                  fac1=10*nroots+1+max(3*(nmax+1)*(mmax+1)*nroots,
     $                                 9+5*nroots)
                  fac2=nroots*3*(jmax+1)*(mmax+1)*((imax+1)*(lmax+1)+
     $                                             (nmax+1))+1
               end if
               fac3=numint+ncontk*npriml
               fac4=numint*2+lenblk*4
c
               space=maxcor-test-1000
               numkl=min(nkl,(space/max(fac1,fac2))/nij)
               if (numkl.le.0) call lnkerr('internal error in driver'//
     $                       ' with numkl')
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
     $                   'error in driver with top1')
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
     $                        (maxcor-h+1)/2)
                  if (lentmp.lt.nijkl) call lnkerr('internal error '//
     $                 'with length of temps for vfmint')
                  temp=h
                  temp1=temp+lentmp
c                 temp=h
c                 temp1=temp+nijkl
c                 tp1=temp1+nijkl
                  top2=max(h+3*(nmax+1)*(jmax+1)*(mmax+1)*lenv,
     $                    temp1+lentmp)
               end if
c
c              if (top2.ne.test+fac2*nijkl) call lnkerr('internal '//
c    $                   ' error in driver with top2')
c
               conint=half+nconti*ncontj*nkl*lenblk
               tmp1=conint+ncint*lenblk
               lent=maxcor-tmp1+1
               if (lent.lt.ncontk*npriml) call lnkerr('internal '//
     $               'error with temporary array for trans')
               top3=tmp1+ncontk*npriml
               if (top3.ne.test+fac3) call lnkerr('internal error'//
     $                   ' in driver with top3')
c
               numint=ncint*lenblk
               plbls=half
               flbls=plbls+numint*nintgr
c              if (flbls.gt.conint) call lnkerr('internal error in '//
c    $                                  'driver with flbls')
               if (flbls+lenblk*4.le.conint) then
                  top4=conint+numint
               else
                  flbls=conint+numint
                  top4=flbls+lenblk*4
               end if
c              if (top4.ne.test+fac4) call lnkerr('internal error'//
c    $                   ' in driver with top4')
c
               top=max(top1,top2,top3,top4)
               maxc=max(top,maxc)
               maxc1=max(maxc1,top1)
               maxc2=max(maxc2,top2)
               maxc3=max(maxc3,top3)
               maxc4=max(maxc4,top4)
               if (top.gt.maxcor) call lnkerr('driver core')
c
c              --- zero half-transformed integrals ---
c               write(iout,*) nconti,ncontj,nkl,lenblk
c               write(iout,*) 'calling rzero' 
               call rzero(z(half),nconti*ncontj*nkl*lenblk)
c
c              --- get primitive information ---
c               write(iout,*) 'calling ldexp first time' 
               call ldexp(z(a),z(ar),ex,nprim,nij,
     $                    ptprim(symcen(1),angmom(1)),
     $                    ptprim(symcen(1),angmom(1))+nprimi-1,
     $                    ptprim(symcen(2),angmom(2)),
     $                    ptprim(symcen(2),angmom(2))+nprimj-1,
     $                    symcen(1),symcen(3),symcen(1),symcen(2),c,
     $                    z(t1),z(xyza),z(xyzam1),z(xyzam3),nat,nprimi)
c
c               write(iout,*) 'calling ldexp second time' 
               call ldexp(z(b),z(br),ex,nprim,nkl,
     $                    ptprim(symcen(3),angmom(3)),
     $                    ptprim(symcen(3),angmom(3))+nprimk-1,
     $                    ptprim(symcen(4),angmom(4)),
     $                    ptprim(symcen(4),angmom(4))+npriml-1,
     $                    symcen(1),symcen(3),symcen(3),symcen(4),c,
     $                    z(t1),z(xyzb),z(xyzbm1),z(xyzbm3),nat,nprimk)
c
c              --- loop over sets of kl primitives ---
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
c                  write(iout,*) 'calling prims' 
                  call prims(z(a),z(ar),z(xyza),z(xyzam1),z(xyzam3),
     $                       nij,z(b),z(br),z(xyzb),z(xyzbm1),
     $                       z(xyzbm3),nkl,minkl,maxkl,z(f00),z(b00),
     $                       z(b10),z(bp01),z(c00),z(cp00),
     $                       z(ab),z(aplusb),z(urho),z(wt),z(denom),
     $                       z(test),z(t1),z(t2),z(t3),z(t4),z(t5),
     $                       z(t6),z(t7),z(t8),z(t9),nijkl,lenv,nmax,
     $                       mmax,nroots,nbatch,cutexp,pi252)
c
                   ntotal=ntotal+nijkl*lenblk
                   nonzer=nonzer+nbatch*lenblk
c
c                  --- if [ss;ss] block, f00, which is overlapped with 
c                      prmint is the primitive, so go to
c                      contraction.
               if (nmax+mmax.eq.0) then
c                  write(iout,*) 'calling trssss' 
                  call trssss(z(half),z(temp1),z(f00),z(test),
     $                        cont(ptcont(symcen(1),angmom(1))),
     $                        cont(ptcont(symcen(2),angmom(2))),
     $                        nijkl,nprimi,nconti,nprimj,ncontj,
     $                        numkl,minkl,nbatch)
               else
c
c              --- form two-dimensional integrals ---
               rbatch=nbatch*nroots
c                  write(iout,*) 'calling vmkghi' 
               call vmkghi(z(g),z(h),z(i2),z(f00),z(b00),z(b10),
     $                     z(bp01),z(c00),z(cp00),nmax,mmax,imax,jmax,
     $                     kmax,lmax,c,nat,symcen,rbatch,lenv)
c
c              --- form primitive integrals ---
c                  write(iout,*) 'calling vfmint' 
               call vfmint(z(i2),z(half),z(temp),lenblk,rbatch,
     $                     angmom,imax,jmax,mmax,lmax,npint,nroots,
     $                     nx,ny,nz,lenxyz,mintyp,nocart,nbtype,minkl,
     $                     nijkl,z(test),nbatch,
     $                     cont(ptcont(symcen(1),angmom(1))),
     $                     cont(ptcont(symcen(2),angmom(2))),
     $                     nprimi,nconti,nprimj,ncontj,nkl,numkl,
     $                     z(temp1),nij,lentmp)
               end if
c
 1001          continue
 1002          continue
c
c              --- contract the primitive integrals ---
c                  write(iout,*) 'calling trans' 
               call trans(z(half),z(conint),nprimk,npriml,
     $                    nconti*ncontj*lenblk,ncontk,ncontl,z(tmp1),
     $                     lent,
     $                     cont(ptcont(symcen(3),angmom(3))),
     $                     cont(ptcont(symcen(4),angmom(4))))
c                  write(iout,*) 'exiting trans'
c
c
c              call trans2(z(prmint),z(half),z(conint),nprimi,
c    $                     nprimj,nprimk,npriml,nconti,ncontj,
c    $                     ncontk,ncontl,
c    $                     cont(ptcont(symcen(1),angmom(1))),
c    $                     cont(ptcont(symcen(2),angmom(2))),
c    $                     cont(ptcont(symcen(3),angmom(3))),
c    $                     cont(ptcont(symcen(4),angmom(4))),
c    $                     z(tmp1),z(tmp2),lenblk)
c
               if (nintgr.eq.1) then
c
c                 --- 64 bit machines ---
                  call out64(z(conint),nconti,ncontj,ncontk,ncontl,
     $                       lenblk,symcen,angmom,start,nat,nbtype,nnp,
     $                       ints,labels,lenbuf,numbuf,
     $                       ncint,z(flbls),numint,z(plbls),
     $                       nbasis,nobf,nactul,drop,pkindx)
               else
c                  write(iout,*) 'calling out32' 
c
c                 --- 32 bit machines ---
                  call out32(z(conint),nconti,ncontj,ncontk,ncontl,
     $                       lenblk,symcen,angmom,start,nat,nbtype,nnp,
     $                       ints,labels,lenbuf,numbuf,
     $                       ncint,z(flbls),numint,z(plbls),
     $                       nbasis,nobf,nactul,drop,pkindx)
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
c     --- flush final integral buffer ---
      nactul=nactul+numbuf-1
      ntotc=nnp*(nnp+1)/2
      ints(1)=-numbuf
      call iosys('write real "unsorted ao integrals" on rints '//
     $     'without rewinding',lenbuf,labels,0,' ')
      call iosys('write real "unsorted ao integrals" on rints '//
     $     'without rewinding',lenbuf,ints,0,' ')
      call iosys('endfile "unsorted ao integrals" on rints',0,0,0,' ')
      if(prnt) then
         pcnt1=float(nonzer)*100.0/float(ntotal)
         pcnt2=float(nactul)*100.0/float(ntotc)
c        write(iout,9010) maxc,maxc1,maxc2,maxc3,maxc4
         write(iout,9020) ntotal,nonzer,pcnt1,ntotc,nactul,pcnt2
      endif
c
c
      call getmem(-ngot,p,idum,'driver',idum)
      return
      end
