*deck @(#)driver.f	5.2  2/5/95
      subroutine driver(ptprim,noprim,nbtype,ex,c,
     #                  nx,ny,nz,
     #                  lenxyz,nobf,nocart,mintyp,minmom,maxmom,
     #                  z,nat,nprim,nderiv,ops,cutexp,acore,
     #                  prnt,nocont,ptcont,
     #                  cont,ncont,start,ints,labels,lenbuf,nbf,nnp,
     #                  cutoff,last,numbuf,lenb2)
c
c***begin prologue     driver
c***date written       891208   (yymmdd)
c
c***revision date      910705   (yymmdd)
c   5 july    1991     rlm at lanl
c      32 bit version.
c***keywords
c***author             saxe, paul (lanl) and lengsfield, byron  (llnl)
c***source             @(#)driver.f	5.2   2/5/95
c***purpose
c         computes derivative integrals over contracted gaussians
c         and then writes a chained file to disk
c***description
c         this is a modified version m313.f which allows subsequent
c         codes to process one set of derivative integrals at a time
c***references         (none)
c
c***routines called    (none)
c
c***end prologue       driver
c
      implicit integer(a-z)
c
      character*(*) ops
      real*8 ex(nprim),c(3,nat),z(*)
      real*8 cutexp,pi252
      real*8 cont(ncont),ints(lenbuf,3*nat),cutoff
      integer last(3*nat),numbuf(3*nat)
      integer acore(*),labels(lenb2,3*nat)
      integer ptprim(nat,nbtype),noprim(nat,nbtype)
      integer nobf(nbtype)
      integer nocart(nbtype),mintyp(nbtype)
      integer minmom(nbtype),maxmom(nbtype)
      integer nx(lenxyz),ny(lenxyz),nz(lenxyz)
      integer symcen(4),angmom(4),dercen(4),dermom(4)
      integer nocont(nat,nbtype),ptcont(nat,nbtype),start(nat,nbtype)
      logical pikcen,pikang,prnt
      logical debug
c
      parameter (debug=.false.)
c
      common /io/     inp,iout
c
c     ----- get 2*pi**(5/2) -----
c
      call iosys('read real pi from rwf',1,pi252,0,' ')
      pi252=2.0d+00*pi252**2.5d+00
c
c     ----- initialize counters -----
c
      maxc=0
      maxc1=0
      maxc2=0
      next=0
      nder=3*nat
      do 111 i=1,nder
         numbuf(i)=1
         last(i)=-1
111   continue
      ntotal=0
      nonzer=0
      nactul=0
c
c
c     call getscm(0,z,canget,'m323:driver',0)
c     call iosys('read integer maxsiz from rwf',1,maxsiz,0,' ')
c
c      if (need1.gt.canget) then
c      write (iout,191) maxsiz,canget,need1,need2
c  191 format(/,' driver: possible field length:',i9,/,
c     #         ' which translates to maxcor of:',i9,
c     #       /,'        minimum core needed is:',i9,/,
c     #         '       maximum core useable is:',i9)
c      call lnkerr('not enough core available for m323')
c      end if
c
c      needed=min(canget,need2)
c
c     ----- get all the memory that is available -----
c
      call getscm(0,z,maxcor,'m323:driver',0)
c
      cencod=-1
c
c     ----- loop over possible centres -----
c
    1 continue
c
         if (pikcen(symcen,nat,cencod)) go to 100
c
c        ----- fix the order of centres in the integrals -----
c
         call redund(symcen,angmom,dercen,dermom,npass)
      if (npass.eq.0) go to 1
c
c     ----- loop over angular momenta of functions -----
c
         angcod=-1
    2    continue
c
      if (pikang(angmom,symcen,noprim,nat,nbtype,cencod,
     #                 angcod)) go to 1
c
c              ----- order centres and number of derivatives -----
c
               call redund(symcen,angmom,dercen,dermom,npass)
c
               imax=maxmom(dermom(1))
               jmax=maxmom(dermom(2))
               kmax=maxmom(dermom(3))
               lmax=maxmom(dermom(4))
               nmax=imax+jmax
               mmax=kmax+lmax
               nroots=(nmax+mmax+nderiv)/2+1
               nprimi=noprim(dercen(1),dermom(1))
               nprimj=noprim(dercen(2),dermom(2))
               nprimk=noprim(dercen(3),dermom(3))
               npriml=noprim(dercen(4),dermom(4))
               nconti=nocont(dercen(1),dermom(1))
               ncontj=nocont(dercen(2),dermom(2))
               ncontk=nocont(dercen(3),dermom(3))
               ncontl=nocont(dercen(4),dermom(4))
               icpt=ptcont(dercen(1),dermom(1))
               jcpt=ptcont(dercen(2),dermom(2))
               kcpt=ptcont(dercen(3),dermom(3))
               lcpt=ptcont(dercen(4),dermom(4))
               nfi=nocart(dermom(1))
               nfj=nocart(dermom(2))
               nfk=nocart(dermom(3))
               nfl=nocart(dermom(4))
               nij=nprimi*nprimj
               nkl=nprimk*npriml
               npint=nij*nkl
               lenblk=nfi*nfj*nfk*nfl
               numint=npint*lenblk
               ntotal=ntotal+numint
               ncint=nconti*ncontj*ncontk*ncontl
c
c              ----- for derivatives, the number of centres to
c                                           differentiate
c
               if (nderiv.eq.0.or.npass.eq.0) then
                  ndcen=0
               else
                  if (npass.eq.1) then
                     ndcen=3
                  else if (npass.eq.4) then
                     ndcen=1
                  else
                     ndcen=2
                  end if
               end if
c
c     ----- core allocation. nb there is a great deal of overlapping --
c
               cint=1
               cdint=cint+ncint*lenblk
               ijindx=wpadti(cdint+ncint*lenblk*3*ndcen)
               klindx=ijindx+nij*2
               aij=iadtwp(klindx+nkl*2)
               ar=aij+nij
               xyza=ar+nij
               xyzam1=xyza+nij*3
               xyzam3=xyzam1+nij*3
               bkl=xyzam3+nij*3
               br=bkl+nkl
               xyzb=br+nkl
               xyzbm1=xyzb+nkl*3
               xyzbm3=xyzbm1+nkl*3
               tmptop=xyzbm3+nkl*3
c
c              ----- work out how long the vectors can be -----
c
               if (nmax+mmax.lt.0) then
                  fac1=10+5*nroots
                  fac2=3+3*nroots
               else
                  if (nderiv.eq.0.or.npass.eq.0) then
                     fac1=6+lenblk+10*nroots+
     #                   max(3*(nmax+1)*(mmax+1)*nroots,
     #                         12+5*nroots)
                     fac2=6+lenblk+
     #                   nroots*3*(jmax+1)*(mmax+1)*((imax+1)*(lmax+1)+
     #                                               (nmax+1))
                  else
                     fac1=6+lenblk+nroots*(10+2*ndcen)+
     #                   max(3*(nmax+1)*(mmax+1)*nroots*(1+ndcen),
     #                         12+5*nroots)
                     fac2=max(6+lenblk+
     #                   nroots*3*(imax+1)*(jmax+1)*(mmax+1)*(lmax+1)*
     #                                             (ndcen+1)+
     #                   nroots*3*(nmax+1)*(jmax+1)*(mmax+1)*(ndcen+1),
     #                   fac1)+nroots*3*ndcen
                  end if
               end if
c
               space=iadtwp(maxcor)-tmptop-100
               len=min(nij*nkl,space/max(fac1,fac2))
c
               if(debug) then
                  write(iout,*) 'driver:,space,maxcor,tmptop',
     $                           space,iadtwp(maxcor),tmptop
                  write(iout,*) 'driver:nij,nkl,fac1,fac2,len',
     $                           nij,nkl,fac1,fac2,len
               endif
               if (len.lt.nij) call lnkerr('internal error in driver'//
     #                       ' with len and nij')
               lenv=len*nroots
c
               index=wpadti(tmptop)
               temp=iadtwp(index+len*6)
               if (nderiv.eq.0.or.npass.eq.0) then
                  dg=temp
                  g=dg
               else
                  dg=temp
                  g=dg+lenv*3*(nmax+1)*(mmax+1)*ndcen
               end if
               ab=dg
               aplusb=ab+len
               urho=aplusb+len
               wt=urho+lenv
               denom=wt+lenv
               expon=denom+lenv
               a=expon+len
               b=a+len
               rho=b+len
               t1=rho+len
               t2=t1+len
               t3=t2+len
               t4=t3+len
               t5=t4+lenv
               t6=t5+lenv
               t7=t6+len
               t8=t7+len
c
               if (nderiv.eq.0.or.npass.eq.0) then
                  rhotsq=1
                  alpha=1
               else
                  rhotsq=t1
                  camcb=rhotsq+len*nroots
                  alpha=camcb+len*3
               end if
c
               if (nmax+mmax.lt.0) then
                  f00=denom
                  b00=1
                  b10=1
                  bp01=1
                  c00=1
                  cp00=1
                  top1=wpadti(t8+len)
               else
                  f00=max(g+3*(nmax+1)*(mmax+1)*lenv,t8+len)
                  b00=f00+lenv
                  b10=b00+lenv
                  bp01=b10+lenv
                  c00=bp01+lenv
                  cp00=c00+lenv*3
                  top1=wpadti(cp00+lenv*3)
                  if (nderiv.ne.0.and.npass.ne.0) then
                     dc00=iadtwp(top1)
                     dcp00=dc00+nroots*len*ndcen
                     top1=wpadti(dcp00+nroots*len*ndcen)
                  end if
               end if
c
cps               if (top1.ne.tmptop+fac1*len) call lnkerr('internal '//
cps     #                   'error in driver with top1')
c
               if (nmax+mmax.lt.0) then
                  temp1=ab
                  i2=1
                  h=1
                  temp=1
                  tp1=1
                  top2=wpadti(f00+lenv)
               else
                  if (nderiv.eq.0.or.npass.eq.0) then
                     i2=g
                     h=i2+3*(imax+1)*(jmax+1)*(mmax+1)*(lmax+1)*lenv
                     top2=wpadti(h+3*(nmax+1)*(jmax+1)*(mmax+1)*lenv)
                     pint=h
                     lnpint=iadtwp(maxcor)-pint
                  else
                     di=dg
                     i2=di+nroots*len*3*(imax+1)*(jmax+1)*(mmax+1)*
     #                                  (lmax+1)*ndcen
                     dh=i2+nroots*len*3*(imax+1)*(jmax+1)*(mmax+1)*
     *                                  (lmax+1)
                     h=dh+nroots*len*3*(nmax+1)*(jmax+1)*(mmax+1)*ndcen
                     d1exp=max(iadtwp(top1),
     #                     h+nroots*len*3*(nmax+1)*(jmax+1)*(mmax+1))
                     top2=wpadti(d1exp+nroots*len*3*ndcen)
                     pint=h
                     lnpint=iadtwp(maxcor)-pint
                  end if
               end if
c
cps               if (top2.ne.tmptop+fac2*len) call lnkerr('internal '//
cps     #                   ' error in driver with top2')
               if (lnpint.lt.nij) call lnkerr('not enough space for '//
     #                                        'primitive integrals')
c
               top=max(top1,top2)
               maxc=max(top,maxc)
               maxc1=max(maxc1,top1)
               maxc2=max(maxc2,top2)
               if(debug) then
                  write(iout,*) 'driver:top1,top2,maxc,maxc1,maxc2,'
     $                  //'maxcor',top1,top2,maxc,maxc1,maxc2,maxcor
               endif
               if (top.gt.maxcor) then
                  write(iout,*) 'driver:top1,top2,maxc,maxc1,maxc2,'
     $                 //'maxcor',top1,top2,maxc,maxc1,maxc2,maxcor
                  call lnkerr('driver core')
               endif
c
c              ----- get primitive information -----
c
               call ldexp(z(aij),z(ar),ex,nprim,nij,
     #                    ptprim(dercen(1),dermom(1)),
     #                    ptprim(dercen(1),dermom(1))+nprimi-1,
     #                    ptprim(dercen(2),dermom(2)),
     #                    ptprim(dercen(2),dermom(2))+nprimj-1,
     #                    dercen(1),dercen(3),dercen(1),dercen(2),c,
     #                    z(t1),z(xyza),z(xyzam1),z(xyzam3),nat,nprimi,
     #                    acore(ijindx))
c
               call ldexp(z(bkl),z(br),ex,nprim,nkl,
     #                    ptprim(dercen(3),dermom(3)),
     #                    ptprim(dercen(3),dermom(3))+nprimk-1,
     #                    ptprim(dercen(4),dermom(4)),
     #                    ptprim(dercen(4),dermom(4))+npriml-1,
     #                    dercen(1),dercen(3),dercen(3),dercen(4),c,
     #                    z(t1),z(xyzb),z(xyzbm1),z(xyzbm3),nat,nprimk,
     #                    acore(klindx))
c
c              ----- loop over sets of kl primitives -----
c
               kl=0
c
c              ----- zero space for contracted integrals -----
c
               call rzero(z(cint),ncint*lenblk)
               call rzero(z(cdint),ncint*lenblk*3*ndcen)
c
 2000       continue
c
c                 ----- form the exponential prefactor and toss out
c                        small integrals
c
                  call prefac(z(ar),z(aij),nij,z(br),z(bkl),nkl,
     #                        z(expon),acore(index),acore(ijindx),
     #                        acore(klindx),kl,len,nv,cutexp,pi252,
     #                        z(t1))
c
c..bhl
            if(nv.eq.0) go to 99
c..bhl
c
c                 ----- form auxiliary arrays from primitive info -----
c
                  call prims(z(aij),z(xyza),z(xyzam1),z(xyzam3),
     #                       nij,z(bkl),z(xyzb),z(xyzbm1),
     #                       z(xyzbm3),nkl,z(f00),z(b00),
     #                       z(b10),z(bp01),z(c00),z(cp00),
     #                       z(ab),z(aplusb),z(urho),z(wt),z(denom),
     #                       z(a),z(b),z(rho),z(expon),
     #                       z(t1),z(t2),z(t3),z(t4),z(t5),
     #                       z(t6),z(t7),z(t8),nv,lenv,nmax,
     #                       mmax,nroots,acore(index),len)
c
c
c                 ----- for derivatives, form an array of exponents ---
c
                  if (nderiv.ne.0.and.npass.ne.0) then
                     call getalp(ex(ptprim(dercen(1),dermom(1))),
     #                           ex(ptprim(dercen(2),dermom(2))),
     #                           ex(ptprim(dercen(3),dermom(3))),
     #                           ex(ptprim(dercen(4),dermom(4))),
     #                           acore(index),z(alpha),nv,len,
     #                           z(camcb),z(xyza),z(xyzb),nij,nkl)
                  end if
                  nonzer=nonzer+nv*lenblk
c
c                 ----- for derivative calculations, form needed
c                        primitive derivative entities

                  if (nderiv.ne.0.and.npass.ne.0) then
                     call dprim(nroots,nv,z(rhotsq),z(urho),z(rho),
     #                          z(dc00),z(dcp00),z(alpha),z(a),z(b),
     #                          z(d1exp),c,z(camcb),nat,dercen,
     #                          ndcen,npass,nmax,mmax)
                  end if
c
c     ----- if [ss;ss] block, f00, which is overlapped with prmint ----
c               is the primitive integrals, so go to contraction.
c
               if (nmax+mmax.lt.0) then
c                 call trssss(z(half),z(temp1),z(f00),z(test),
c    #                        cont(ptcont(dercen(1),dermom(1))),
c    #                        cont(ptcont(dercen(2),dermom(2))),
c    #                        nijkl,nprimi,nconti,nprimj,ncontj,
c    #                        numkl,minkl)
               else
               if (nderiv.eq.0.or.npass.eq.0) then
c
c                 ----- form two-dimensional integrals -----
c
                  call vmkghi(z(g),z(h),z(i2),z(f00),z(b00),z(b10),
     #                     z(bp01),z(c00),z(cp00),nmax,mmax,imax,jmax,
     #                     kmax,lmax,c,nat,dercen,nv*nroots)
               else
c
c                 ----- form derivative two-dimensional integrals -----
c
                  call mkdghi(z(g),z(h),z(i2),z(f00),z(b00),z(b10),
     #                        z(bp01),z(c00),z(cp00),nmax,mmax,
     #                        imax,jmax,kmax,lmax,c,nat,dercen,
     #                        nv*nroots,z(dg),z(dh),z(di),
     #                        z(dc00),z(dcp00),ndcen,z(d1exp),npass)
c
c                 ----- form energy and derivative contribution -----
c
                  leni=nroots*nv*(imax+1)*(jmax+1)*(mmax+1)*(lmax+1)
                  call fmdint(z(i2),lenblk,nv*nroots,
     #                        dermom,imax,jmax,mmax,lmax,nroots,
     #                        nx,ny,nz,lenxyz,mintyp,nocart,nbtype,
     #                        nv,z(di),ndcen,leni,
     #                        z(pint),lnpint,cont(icpt),
     #                        nprimi,nconti,cont(jcpt),nprimj,ncontj,
     #                        cont(kcpt),nprimk,ncontk,cont(lcpt),
     #                        npriml,ncontl,acore(index),len,z(cdint),
     #                        ncint)
               end if
               end if
c..bhl
  99       continue
c..bhl
c
            if (kl.lt.nkl) go to 2000
c
c           ----- fix up derivatives of [ii;jj] and [ij;ij] -----
c
c..bhl
            if(ncint*lenblk.ne.0) then
c..bhl
               call fixup(npass,dercen,z(cdint),ncint*lenblk,ndcen)
c
c              ----- transfer the contracted integrals to buffers -----
c
               flbls=wpadti(cdint+ncint*lenblk*3*ndcen)
               plbls=flbls+lenblk*4
c
c              ---- on 32/64 bit machines, 'plbls' is 2 x * array
c
               if (wptoin(1).eq.2) then
                  top5=plbls+2*ncint*lenblk
               else
                  top5=plbls+ncint*lenblk
               end if
c
               if (top5.gt.maxcor) call lnkerr(' plbls')
c
               if (wptoin(1).eq.2) then
                  call dout32(z(cint),nconti,ncontj,ncontk,ncontl,
     #                    lenblk,dercen,dermom,start,nat,nbtype,nnp,
     #                    ints,labels,lenbuf,numbuf,
     #                    ncint,acore(flbls),ncint*lenblk,acore(plbls),
     #                    nbf,nobf,nactul,z(cdint),ndcen,
     #                    npass,last,nder,next)
               else
                  call dout64(z(cint),nconti,ncontj,ncontk,ncontl,
     #                    lenblk,dercen,dermom,start,nat,nbtype,nnp,
     #                    ints,labels,lenbuf,numbuf,
     #                    ncint,acore(flbls),ncint*lenblk,acore(plbls),
     #                    nbf,nobf,nactul,z(cdint),ndcen,
     #                    npass,last,nder,next)
               endif
c..bhl
           end if
c..bhl
c
c     ----- loop back up for the next angular momentum -----
c
            go to 2
c
  100 continue
c
c     ----- flush the final integral buffer -----
c
      do 101 i=1,nder
         nactul=nactul+numbuf(i)-1
         ints(1,i)=numbuf(i)
         labels(1,i)=last(i)
         last(i)=next
         next=next+2*lenbuf
         call iosys('write real "unsorted derivative integrals"'
     #            //' on rdints without rewinding',
     $              lenbuf,labels(1,i),0,' ')
         call iosys('write real "unsorted derivative integrals"'
     #            //' on rdints without rewinding',
     $              lenbuf,ints(1,i),0,' ')
 101  continue
      call iosys('endfile "unsorted derivative integrals" on rdints',
     $           0,0,0,' ')
c
      call iosys('write integer nder to rdints',1,nder,0,' ')
      call iosys('write integer last to rdints',nder,last,0,' ')
c
      if (prnt) write (iout,192) nactul
  192 format (5x,'# integrals kept    ',19x,i10)
c
      if(prnt) write (iout,93) ntotal,nonzer,
     $                         float(nonzer)*100.0/float(ntotal)
   93 format(5x,'# integrals possible',19x,i10,
     $      /5x,'# integrals computed',19x,i10,'(',f5.1,'%)')
c
c
      return
      end
