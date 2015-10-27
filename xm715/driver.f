*deck %W%  %G%
      subroutine driver(ptprim,noprim,nbtype,ex,c,
     $     nx,ny,nz,lenxyz,nobf,nocart,mintyp,minmom,maxmom,
     $     z,maxcor,nat,npf,nnprim,nprim,ops,cutexp,acore,
     $     grad,nderiv,prnt,ndmat,alpha,beta,
     $     nshell,pstart,dpr,d2e,nd2e,
     $     nbf,
     $     cont,ncont,ptcont,nocont,start)
c***begin prologue     driver.f
c***date written       851113   
c***revision date      driver.f 
c
c  01 february  1988    bhl at brl
c     c..bhl.nv  denotes changes to skip subroutine calls
c     if no integrals (nv=0) past the prefactor screening
c
c  16 december  1987    bhl at brl
c     core allocation for tpaadm tpacdm & tpccdm
c
c  08 december  1987    bhl at brl
c     splitting fd2int into two routines (fd2ant,fd2bnt)
c     so cos can compile the routines
c
c  20 november 1987    pws at lanl
c     adding arguments from 'bftgrp' to 'start' for reading in the
c     ci or mcscf two-particle density matrices.
c
c  11 november 1987    pws at lanl
c     adding arguments to call to 'dprim' in preparation for second
c     derivatives.
c
c   3 july  1987    pws at lanl
c     switching density matrix work to in-core to save on ridiculous
c     i/o charges accumulating.
c
c  22 march 1987    pws at lanl
c      creating a 'sizer' routine which figures out how much core
c      we'll need.
c
c
c***keywords
c***author             saxe, paul    (lanl)
c***source             %W%   %G%
c***purpose
c***description
c***references         (none)
c
c***routines called    (none)
c
c***end prologue       driver.f
      implicit none
c     --- input variables -----
      character*(*) ops
      logical prnt
      integer nbtype,lenxyz,maxcor,nat,npf,nnprim,nprim
      integer nderiv,ndmat,nshell,nd2e,nbf
      integer ncont
      real*8 cutexp
c     --- input arrays (unmodified) ---
      integer ptcont(nat,nbtype)
      integer nocont(nat,nbtype)
      integer start(nat,nbtype)
      integer pstart(nat,nbtype)
      integer ptprim(nat,nbtype),noprim(nat,nbtype)
      integer nobf(nbtype)
      integer nocart(nbtype),mintyp(nbtype)
      integer minmom(nbtype),maxmom(nbtype)
      integer nx(lenxyz),ny(lenxyz),nz(lenxyz)
      real*8 cont(ncont)
      real*8 ex(nprim),c(3,nat)
      real*8 alpha(nshell,nshell),beta(nshell,nshell)
      real*8 dpr(nnprim,ndmat)
c     --- input arrays (scratch) ---
      integer acore(*)
      real*8 z(maxcor)
c     --- output arrays ---
      real*8 grad(3,nat)
      real*8 d2e(nd2e)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer symcen(4),angmom(4),dercen(4),dermom(4)
      integer maxc,maxc1,maxc2,numbuf,ntotal,nonzer,nactul
      integer igrp,iatom,itype,jgrp,jatom,jtype,jtmax,ijgrp
      integer kgrp,katom,ktype,ktmax,klgrp
      integer latmax,ltmax,lgrp,latom,ltype
      integer npass,flip
      integer imax,jmax,kmax,lmax,nmax,mmax,nroots
      integer nprimi,nprimj,nprimk,npriml,nconti,ncontj,ncontk,ncontl
      integer icpt,jcpt,kcpt,lcpt,ipstrt,jpstrt,kpstrt,lpstrt
      integer istart,jstart,kstart,lstart
      integer nfi,nfj,nfk,nfl,nij,nkl,npint,lenblk,numint
      integer ndcen,nd2
      integer dij,dkl,ijindx,klindx,aij,ar,xyza,xyzam1,xyzam3
      integer bkl,br,xyzb,xyzbm1,xyzbm3,tmptop
      integer fac1,fac2
      integer space,facmax,lenmin,mincor,len,lenv,index,twopdm
      integer g,dg,d2g
      integer ab,aplusb,urho,wt,denom,expon,a,b,rho,t1,t2,t3,t4,t5
      integer t6,t7,t8
      integer rhotsq,expnts,camcb
      integer f00,b00,b10,bp01,c00,cp00,top1
      integer dc00,dcp00,temp1,i2,h,temp,tp1,top2
      integer di,d1exp,d2exp
      integer dh,d2i,d2h
      integer kl,minkl,nv,leni
      integer iadtwp,wpadti,wptoin
      integer maxtop,top
      real*8 etot,energy,enaa,enac,encc,del
      real*8 pi252
      real*8 ld2e(78)
      real*8 der(3,4)
      logical debug
c
      real*8 e1,e2,e1702
      real*8 small
      logical logkey
c
c
      integer mynodeid,nodeid,nprocs,nnodes,mdtob,mitob,nxtask,next,
     $     ijblk
      include 'msgtypesf.h'
      logical ispar
      common /tcgmesa/ ispar
c
c
      parameter (debug=.false.)
      common /io/     inp,iout
c
c
      mynodeid=nodeid()
      nprocs=nnodes()
c
c     --- get the primitive density matrices for hf or mcscf ---
      if (mynodeid .eq. 0) then
         call iosys('read real "hf primitive density" from rwf',
     $        nnprim*ndmat,dpr,0,' ')
      endif
      if (ispar) call brdcst(101+MSGDBL,dpr,mdtob(nnprim*ndmat),0)
c
c     --- get 2*pi**(5/2) ---
      if (mynodeid .eq.0) then
         call iosys('read real pi from rwf',1,pi252,0,' ')
      endif
      if (ispar) call brdcst(102+MSGDBL,pi252,mdtob(1),0)
      pi252=2.0d+00*pi252**2.5d+00
c
c     --- zero energy, gradients, etc. ---
      etot=0.d0
      energy=0.0d+00
      if (nderiv.ge.1) call rzero(grad,3*nat)
      if (nderiv.ge.2) call rzero(d2e,nd2e)
c
c     --- initialize counters ---
      maxc=0
      maxc1=0
      maxc2=0
      numbuf=1
      ntotal=0
      nonzer=0
      nactul=0
c
c     --- check amount of core needed ---
      call sizer(noprim,nbtype,
     $           nx,ny,nz,
     $           lenxyz,nobf,nocart,mintyp,minmom,maxmom,
     $           nat,npf,nnprim,nprim,ops,cutexp,
     $           nderiv,prnt,ndmat,
     $           nshell,mincor)
c
c     --- let's go ---
      igrp=0
      ijblk=0
      next=nxtask(nprocs)+1
      do 8000 iatom=1,nat
         symcen(1)=iatom
         do 7000 itype=1,nbtype
            if (noprim(iatom,itype).le.0) go to 7000
            angmom(1)=itype
            igrp=igrp+1
            jgrp=0
            do 6000 jatom=1,iatom
               symcen(2)=jatom
               if (jatom.eq.iatom) then
                  jtmax=itype
               else
                  jtmax=nbtype
               end if
               do 5000 jtype=1,jtmax
                  if (noprim(jatom,jtype).le.0) go to 5000
                  angmom(2)=jtype
                  jgrp=jgrp+1
                  ijgrp=igrp*(igrp-1)/2+jgrp
                  ijblk=ijblk+1
                  if (ijblk .ne. next) goto 5000
                  kgrp=0
                  do 4000 katom=1,iatom
                     symcen(3)=katom
                     if (katom.eq.iatom) then
                        ktmax=itype
                     else
                        ktmax=nbtype
                     end if
                     do 3000 ktype=1,ktmax
                        if (noprim(katom,ktype).le.0) go to 3000
                        angmom(3)=ktype
                        if (katom.eq.iatom.and.ktype.eq.itype) then
                           latmax=jatom
                        else
                           latmax=katom
                        end if
                        kgrp=kgrp+1
                        lgrp=0
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
                              if (noprim(latom,ltype).le.0) go to 1000
                              angmom(4)=ltype
                              lgrp=lgrp+1
                              klgrp=kgrp*(kgrp-1)/2+lgrp
c
               if(debug) then 
                  write(iout,77001) iatom,jatom,katom,latom,
     $                 angmom(1),angmom(2),angmom(3),angmom(4),
     $                 igrp,jgrp,kgrp,lgrp 
77001 format(/,' iatm jatm katm latm iang   jang   kang   lang   ',
     $       8i8,/,' igrp jgrp kgrp lgrp ',8i8)
              endif
c
c              --- order centres and number of derivatives ---
               call redund(symcen,angmom,dercen,dermom,npass,flip)
c
c              --- set up variables for this shell block combination.
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
               ipstrt=pstart(dercen(1),dermom(1))
               jpstrt=pstart(dercen(2),dermom(2))
               kpstrt=pstart(dercen(3),dermom(3))
               lpstrt=pstart(dercen(4),dermom(4))
               istart=start(dercen(1),dermom(1))
               jstart=start(dercen(2),dermom(2))
               kstart=start(dercen(3),dermom(3))
               lstart=start(dercen(4),dermom(4))
               nfi=nocart(dermom(1))
               nfj=nocart(dermom(2))
               nfk=nocart(dermom(3))
               nfl=nocart(dermom(4))
               nij=nprimi*nprimj
               nkl=nprimk*npriml
               npint=nij*nkl
               lenblk=nfi*nfj*nfk*nfl
               numint=npint*lenblk
c
c              --- for derivatives, the number of centres to
c                  differentiate
               if (nderiv.eq.0.or.npass.eq.0) then
                  ndcen=0
                  nd2=0
               else
                  if (npass.eq.1) then
                     ndcen=3
                     nd2=6
                  else if (npass.eq.4) then
                     ndcen=1
                     nd2=1
                  else
                     ndcen=2
                     nd2=3
                  end if
               end if
               ntotal=ntotal+numint
c
c              --- core allocation. nb there is a great deal of overlapping 
               dij=1
               dkl=dij+nprimi*nprimj*nfi*nfj*ndmat
               ijindx=wpadti(dkl+nprimk*npriml*nfk*nfl*ndmat)
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
               tmptop=wpadti(xyzbm3+nkl*3)
c
c              --- work out how long the vectors can be ---
               if (nmax+mmax.lt.0) then
                  fac1=10+5*nroots
                  fac2=3+3*nroots
               else
                  if (nderiv.eq.0.or.npass.eq.0) then
                     fac1=6+wptoin(lenblk+10*nroots+
     $                   max(3*(nmax+1)*(mmax+1)*nroots,
     $                         12+5*nroots))
                     fac2=6+wptoin(lenblk+
     $                   nroots*3*(jmax+1)*(mmax+1)*((imax+1)*(lmax+1)+
     $                                               (nmax+1)))
                  else if (nderiv.eq.1) then
                     fac1=6+wptoin(lenblk+nroots*(10+2*ndcen)+
     $                   max(3*(nmax+1)*(mmax+1)*nroots*(1+ndcen),
     $                         12+5*nroots))
                     fac2=max(6+wptoin(lenblk+
     $                   nroots*3*(imax+1)*(jmax+1)*(mmax+1)*(lmax+1)*
     $                                             (ndcen+1)+
     $                   nroots*3*(nmax+1)*(jmax+1)*(mmax+1)*(ndcen+1)),
     $                   fac1)+wptoin(nroots*3*ndcen)
                  else if (nderiv.eq.2) then
                     fac1=6+wptoin(lenblk+nroots*(10+2*ndcen)+
     $                   max(3*(nmax+1)*(mmax+1)*nroots*(1+ndcen+nd2),
     $                         12+5*nroots))
                     fac2=max(6+wptoin(lenblk+
     $                   nroots*3*(imax+1)*(jmax+1)*(mmax+1)*(lmax+1)*
     $                                             (nd2+ndcen+1)+
     $                   nroots*3*(nmax+1)*(jmax+1)*(mmax+1)*
     $                    (nd2+ndcen+1)),fac1)+
     $                    wptoin(nroots*(3*ndcen+nd2))
                  else
                     call plnkerr('too many derivatives asked for',100)
                  end if
               end if
c
c
               space=maxcor-tmptop-100
               facmax=max(fac1,fac2)
               lenmin=facmax*nij
               if(space.lt.lenmin) then
                  mincor=lenmin+tmptop+101
                  if(mincor.gt.maxtop) then
                     call plnkerr('insufficient memory for m715',101)
c                     write(iout,*)' insufficient memory for m715: '
c                     write(iout,*)' ihave  ineed ',maxtop,mincor
c                     write(iout,22011) nij,nkl,fac1,fac2,npint,lenblk
c22011                format(' nij nkl fac1 fac2 npint lenblk',
c     $                      6(2x,i8))
                  end if
                  call getscm(mincor,z,maxcor,' m715:bhl ',0)
                  space=maxcor-tmptop-100
               endif
c
c
               len=min(nij*nkl,space/max(fac1,fac2))
c
               if (len.lt.nij) then
c                  write(iout,*)' len nij nkl ',len,nij,nkl
c                  write(iout,*)' lenmin  facmax ',
c     $                        lenmin,facmax
c                  write(iout,*)' fac1 fac2 space ',fac1,fac2,space
                  call plnkerr('internal error in driver'//
     $                        ' with len and nij',102)
               end if
               lenv=len*nroots
c
               index=tmptop
c
               twopdm=iadtwp(index+len*6)
               if (nderiv.eq.0.or.npass.eq.0) then
                  dg=twopdm+len*lenblk
                  d2g=dg
                  g=d2g
               else if (nderiv.eq.1) then
                  dg=twopdm+len*lenblk
                  d2g=dg
                  g=dg+lenv*3*(nmax+1)*(mmax+1)*ndcen
               else if (nderiv.eq.2) then
                  dg=twopdm+len*lenblk
                  d2g=dg+lenv*3*(nmax+1)*(mmax+1)*ndcen
                  g=d2g+lenv*3*(nmax+1)*(mmax+1)*nd2
               end if
c
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
                  expnts=1
               else
                  rhotsq=t1
                  camcb=rhotsq+len*nroots
                  expnts=camcb+len*3
               end if
c
               if (nmax+mmax.lt.0) then
                  f00=denom
                  b00=1
                  b10=1
                  bp01=1
                  c00=1
                  cp00=1
                  top1=t8+len
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
cps               if (abs(top1-(tmptop+fac1*len)).gt.4) then
cps                  call lnkerr('internal error in driver with top1')
cps               end if
c
               if (nmax+mmax.lt.0) then
                  temp1=ab
                  i2=1
                  h=1
                  temp=1
                  tp1=1
                  top2=f00+lenv
               else
                  if (nderiv.eq.0.or.npass.eq.0) then
                     i2=g
                     h=i2+3*(imax+1)*(jmax+1)*(mmax+1)*(lmax+1)*lenv
                     top2=wpadti(h+3*(nmax+1)*(jmax+1)*(mmax+1)*lenv)
                  else if (nderiv.eq.1) then
                     di=dg
                     i2=di+nroots*len*3*(imax+1)*(jmax+1)*(mmax+1)*
     $                                  (lmax+1)*ndcen
                     dh=i2+nroots*len*3*(imax+1)*(jmax+1)*(mmax+1)*
     *                                  (lmax+1)
                     h=dh+nroots*len*3*(nmax+1)*(jmax+1)*(mmax+1)*ndcen
                     d1exp=max(iadtwp(top1),
     $                     h+nroots*len*3*(nmax+1)*(jmax+1)*(mmax+1))
                     d2exp=d1exp
                     top2=wpadti(d1exp+nroots*len*3*ndcen)
                  else if (nderiv.eq.2) then
                     di=dg
                     d2i=di+nroots*len*3*(imax+1)*(jmax+1)*(mmax+1)*
     $                    (lmax+1)*ndcen
                     i2=d2i+nroots*len*3*(imax+1)*(jmax+1)*(mmax+1)*
     $                    (lmax+1)*nd2
                     d2h=i2+nroots*len*3*(imax+1)*(jmax+1)*(mmax+1)*
     $                    (lmax+1)
                     dh=d2h+nroots*len*3*(nmax+1)*(jmax+1)*(mmax+1)*nd2
                     h=dh+nroots*len*3*(nmax+1)*(jmax+1)*(mmax+1)*ndcen
                     d1exp=max(iadtwp(top1),
     $                    h+nroots*len*3*(nmax+1)*(jmax+1)*(mmax+1))
                     d2exp=d1exp+nroots*len*3*ndcen
                     top2=wpadti(d2exp+nroots*len*nd2)
                  end if
               end if
c
cps               if (abs(top2-(tmptop+fac2*len)).gt.4) then
cps                  call lnkerr('internal error in driver with top2')
cps               end if
c
               top=max(top1,top2)
               maxc=max(top,maxc)
               maxc1=max(maxc1,top1)
               maxc2=max(maxc2,top2)
               if (top.gt.maxcor) then
                 call getscm(top+1000,z,maxcor,'m715:driver',0)
c                 write(iout,*)' m715:  top maxcor ',top,maxcor
c..bhl                call lnkerr('driver core')
               endif
c
c              --- pick up and transform the ao two-particle density
c                  matrix to the primitive basis. first get local
c                  blocks of the density matrix.
               call get1dm(dpr,nnprim,ndmat,z(dij),nprimi,nprimj,
     $                     nfi,nfj,ipstrt,jpstrt)
               call get1dm(dpr,nnprim,ndmat,z(dkl),nprimk,npriml,
     $                     nfk,nfl,kpstrt,lpstrt)
c
c               if (logkey(ops,'m715=print=local-density-matrices',
c     $                   .false.,' ')) then
c                  call prntld(symcen,angmom,nprimi,nprimj,nprimk,
c     $                        npriml,nfi,nfj,nfk,nfl,ipstrt,jpstrt,
c     $                        kpstrt,lpstrt,ndmat,z(dij),z(dkl))
c               end if
c
c              --- get primitive information ---
               call ldexp(z(aij),z(ar),ex,nprim,nij,
     $                    ptprim(dercen(1),dermom(1)),
     $                    ptprim(dercen(1),dermom(1))+nprimi-1,
     $                    ptprim(dercen(2),dermom(2)),
     $                    ptprim(dercen(2),dermom(2))+nprimj-1,
     $                    dercen(1),dercen(3),dercen(1),dercen(2),c,
     $                    z(t1),z(xyza),z(xyzam1),z(xyzam3),nat,
     $                    nprimi,acore(ijindx))
c
               call ldexp(z(bkl),z(br),ex,nprim,nkl,
     $                    ptprim(dercen(3),dermom(3)),
     $                    ptprim(dercen(3),dermom(3))+nprimk-1,
     $                    ptprim(dercen(4),dermom(4)),
     $                    ptprim(dercen(4),dermom(4))+npriml-1,
     $                    dercen(1),dercen(3),dercen(3),dercen(4),c,
     $                    z(t1),z(xyzb),z(xyzbm1),z(xyzbm3),nat,nprimk,
     $                    acore(klindx))
c
c              --- loop over sets of kl primitives ---
               kl=0
               del=0.0d+00
               enaa=0.d0
               enac=0.d0
               encc=0.d0
               call rzero(der,12)
               if (nderiv.ge.2) call rzero(ld2e,78)
c
 200           continue
                  minkl=kl+1
c
c                 --- form the exponential prefactor and toss out
c                     small integrals
                  call prefac(z(ar),z(aij),nij,z(br),z(bkl),nkl,
     $                        z(expon),acore(index),acore(ijindx),
     $                        acore(klindx),kl,len,nv,cutexp,pi252,
     $                        z(t1))
                  if(nv.eq.0) go to 99
c
c                    --- form the two-particle density matrix from the
c                        one-particle
c
                     call gj2pdm(z(twopdm),
     $                    z(dij),z(dkl),
     $                    nprimi,nprimj,nprimk,npriml,
     $                    nfi,nfj,nfk,nfl,ipstrt,jpstrt,kpstrt,
     $                    lpstrt,len*lenblk,acore(index),nv,len,
     $                    alpha,nshell,ndmat)
c
c
c                    --- form auxiliary arrays from primitive info ---
                     call prims(z(aij),z(xyza),z(xyzam1),z(xyzam3),
     $                          nij,z(bkl),z(xyzb),z(xyzbm1),
     $                          z(xyzbm3),nkl,z(f00),z(b00),
     $                          z(b10),z(bp01),z(c00),z(cp00),
     $                          z(ab),z(aplusb),z(urho),z(wt),z(denom),
     $                          z(a),z(b),z(rho),z(expon),
     $                          z(t1),z(t2),z(t3),z(t4),z(t5),
     $                          z(t6),z(t7),z(t8),nv,lenv,nmax,
     $                          mmax,nroots,acore(index),len)
c
c                    --- for derivatives, form an array of exponents ---
                     if (nderiv.ne.0.and.npass.ne.0) then
                        call getalp(ex(ptprim(dercen(1),dermom(1))),
     $                              ex(ptprim(dercen(2),dermom(2))),
     $                              ex(ptprim(dercen(3),dermom(3))),
     $                              ex(ptprim(dercen(4),dermom(4))),
     $                              acore(index),z(expnts),nv,len,
     $                              z(camcb),z(xyza),z(xyzb),nij,nkl)
                     end if
                     if (nderiv.eq.0.or.npass.eq.0) then
                        nonzer=nonzer+nv*lenblk
                     else if (nderiv.eq.1) then
                        nonzer=nonzer+nv*lenblk*(3*ndcen+1)
                     else if (nderiv.eq.2) then
                        if (npass.eq.1) then
                           nonzer=nonzer+nv*lenblk*(45+9+1)
                        else if (npass.eq.2.or.npass.eq.3) then
                           nonzer=nonzer+nv*lenblk*(21+6+1)
                        else if (npass.eq.4) then
                           nonzer=nonzer+nv*lenblk*(6+3+1)
                        else
                           nonzer=nonzer+nv*lenblk
                        end if
                     end if
c
c                    --- for derivative calculations, form needed
c                        primitive derivative entities
c
                     if (nderiv.ne.0.and.npass.ne.0) then
                        call dprim(nroots,nv,z(rhotsq),z(urho),z(rho),
     $                             z(dc00),z(dcp00),z(expnts),z(a),z(b),
     $                             z(d1exp),c,z(camcb),nat,dercen,
     $                             ndcen,npass,nmax,mmax,z(d2exp),nd2,
     $                             nderiv)

                     end if
                     if (nderiv.eq.0.or.npass.eq.0) then
c
c                       --- form two-dimensional integrals ---
                        call vmkghi(z(g),z(h),z(i2),z(f00),z(b00),
     $                              z(b10),z(bp01),z(c00),z(cp00),
     $                              nmax,mmax,imax,jmax,kmax,lmax,
     $                              c,nat,dercen,nv*nroots)
c
c                       --- form energy contribution ---
                        call vfmint(z(i2),z(twopdm),
     $                              lenblk,nv*nroots,
     $                              dermom,imax,jmax,mmax,lmax,nroots,
     $                              nx,ny,nz,lenxyz,mintyp,nocart,
     $                              nbtype,nv,del)
c
                     else if (nderiv.eq.1) then
c
c                       --- form derivative two-dimensional integrals
                        call mkdghi(z(g),z(h),z(i2),z(f00),z(b00),
     $                              z(b10),z(bp01),z(c00),z(cp00),
     $                              nmax,mmax,imax,jmax,kmax,lmax,
     $                              c,nat,dercen,nv*nroots,z(dg),z(dh),
     $                              z(di),z(dc00),z(dcp00),ndcen,
     $                              z(d1exp),npass)
c
c                       --- form energy and derivative contribution ---
                        leni=nroots*nv*3*(imax+1)*(jmax+1)*(mmax+1)*
     $                       (lmax+1)
c
                        call fmdint(z(i2),z(twopdm),
     $                              lenblk,nv*nroots,
     $                              dermom,imax,jmax,mmax,lmax,nroots,
     $                              nx,ny,nz,lenxyz,mintyp,nocart,
     $                              nbtype,nv,del,z(di),ndcen,der,
     $                              leni)
c
                     else if (nderiv.eq.2) then
c
c                       --- form derivative two-dimensional integrals
                        call md2ghi(z(g),z(h),z(i2),z(f00),z(b00),
     $                              z(b10),z(bp01),z(c00),z(cp00),
     $                              nmax,mmax,imax,jmax,kmax,lmax,
     $                              c,nat,dercen,nv*nroots,z(dg),
     $                              z(dh),z(di),
     $                              z(dc00),z(dcp00),ndcen,nd2,z(d1exp),
     $                              z(d2exp),npass,z(d2g),z(d2h),z(d2i))
c
c                       --- form energy and derivative contribution -
                        leni=nroots*nv*3*(imax+1)*(jmax+1)*(mmax+1)*
     $                       (lmax+1)
                        if(nroots.gt.5) then
                           call fd2ant(z(i2),z(twopdm),lenblk,nv*nroots,
     $                                 dermom,imax,jmax,mmax,lmax,
     $                                 nroots,nx,ny,nz,lenxyz,mintyp,
     $                                 nocart,nbtype,nv,del,z(di),ndcen,
     $                                 der,leni,z(d2i),nd2e,ld2e,npass)
                        else
                           call fd2bnt(z(i2),z(twopdm),lenblk,nv*nroots,
     $                                 dermom,imax,jmax,mmax,lmax,
     $                                 nroots,nx,ny,nz,lenxyz,mintyp,
     $                                 nocart,nbtype,nv,del,z(di),ndcen,
     $                                 der,leni,z(d2i),nd2e,ld2e,npass)
                        endif
                     end if
  99              continue
c
               if (kl.lt.nkl) go to 200
c
c              ---
               energy=energy+del
c
               if(debug) then
                  write(iout,11022)nroots,energy,del
     $                         
11022            format(/,' nroots  = ',i5,
     $                  /,' energy  del    ',2(2x,f12.8))

                  small=1.d-08
c                  if(abs(del-(encc+enac+enaa)).gt.small) 
c     $               call lnkerr(' bug ')
               end if
c
c              --- move the momentum group-block gradients to the
c                  total array ---
               if (nderiv.ge.1.and.npass.ne.0) then
                  call movder(der,grad,nat,dercen,npass)
               end if
               if (nderiv.ge.2.and.npass.ne.0) then
                  call movd2e(ld2e,d2e,nd2e,dercen,ndcen,npass)
               end if
c
c
 1000                      continue
 2000                   continue
 3000                continue
 4000             continue
                  next=nxtask(nprocs)+1
 5000          continue
 6000       continue
 7000    continue
 8000 continue
c----
cTCGMSG     BARRIER
c----
      next=nxtask(-nprocs)
c
c
      if (ispar) call dgop(101+MSGDBL,energy,1,'+')
      if (ispar) call dgop(102+MSGDBL,grad,3*nat,'+')
      if (nderiv .ge. 2 .and. ispar) call dgop(103+MSGDBL,d2e,nd2e,'+')
      if (mynodeid .eq. 0) then
         if(prnt) write (iout,93) ntotal,nonzer,
     $        float(nonzer)*100.0/float(ntotal)
 93      format(5x,'#integrals possible',19x,i10,
     $        /5x,'#integrals computed',19x,i10,'(',f5.1,'%)')
c     
c     --- check against dft coulomb energy.
         call iosys('read real "hf 1e energy" from rwf',1,e1,0,' ')
         call iosys('read real "hf coulomb energy" from rwf',
     $        1,e2,0,' ')
         call iosys('read real "hf m702 1e energy" from rwf',
     $        1,e1702,0,' ')

         if(abs((e2-energy)/e2).gt.1.0d-06.and.
     $        abs(energy).gt.1.0d-06) then
            write (iout,87) e1702,e1,energy,e2
 87         format (/5x,'calculated one-electron energy is:',g20.12,
     $           /,5x,'  one-electron energy from scf is:',g20.12,
     $           /,5x,'     calculated coulomb energy is:',g20.12,
     $           /,5x,'       coulomb energy from scf is:',g20.12)
c            call plnkerr('energies do not agree !!!',103)
         end if
      endif
c
c
      return
      end
