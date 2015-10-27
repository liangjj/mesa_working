      subroutine mcdnj(xj,hess,den,aii,aij,
     $     dab,ldab,deg,fab,lfab,noba,lnoba,nco,nao,nvirt,nob,c,isymc,
     $     nsym,noci,temp,tv,itran,thrsh,mmr,mms,ndab,nfiv,nfav,
     $     lhaa,lhac,lhca,lhcc,
     $     xjbuf,nocc,nbf,nbufj,jstrt,jnij,jmij,jnkl,jmkl,jkcore,
     $     bufix,bufax,lbufso)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
c     implicit real*8(a-h,o-p,r-z),             integer*2(q)
c
      dimension xj(2),hess(2),den(2),dab(2),ldab(2),deg(2),fab(2),
     $     lfab(2),nco(2),nao(2),nvirt(2),nob(2),c(2),isymc(2),temp(2),
     $     tv(2),noba(2),lnoba(2),aij(10,2),aii(2),xjbuf(2)
      real*8 bufix(lbufso),bufax(lbufso)
c
      character*3 ians
      dimension cz(2),ihd(14)
c
c
      common /io/ inp,iout
c
      common / number / zero,pt5,one,two,four,eight
c
c----------------------------------------------------
c
c  this program processes the coulomb integral blocks
c  these blocks contribute to the orbital hessian and
c  to the  fab  fock operator. in addition, some of these
c  blocks are needed to compute derivatives of fock matrix
c  elements. these last terms are needed when ci-coupling
c  is invoked and are called i-vector terms in our convention
c---------------------------------------------------------------
c
c
      ikt=0
c
c      do 1000 irskl=1,nklv
c
c      mk     = mklrs( 1,irskl) + 1
c      ml     = mklrs( 2,irskl) + 1
c      mr     = mklrs( 3,irskl) + 1
c      ms     = mklrs( 4,irskl) + 1
c      idtype = mklrs( 5,irskl)
c      nu     = mklrs( 6,irskl)
c      mll    = mklrs( 7,irskl)
c      ipqrs  = mklrs( 8,irskl)
c      itype  = mklrs( 9,irskl)
c      nok    = mklrs(10,irskl)
c      nol    = mklrs(11,irskl)
c      nbr    = mklrs(12,irskl)
c      nbs    = mklrs(13,irskl)
c      mkl    = mklrs(14,irskl)
c      mrs    = mklrs(15,irskl)
c      nrs    = mklrs(16,irskl)
c      nkl    = mklrs(17,irskl)
c      idad0  = mklrs(18,irskl)
c      idad1  = mklrs(19,irskl)
c      idad2  = mklrs(20,irskl)
c
      nnp=nbf*(nbf+1)/2
      nnb=nocc*(nocc+1)/2
      nnbf=nbf*nbf
c
      mk     = 1
      ml     = 1
      mr     = 1
      ms     = 1
      idtype = 1
      nu     = 0
      mll    = 0
      ipqrs  = 1
      itype  = 1
      nok    = nocc
      nol    = nocc
      nbr    = nbf
      nbs    = nbf
      mkl    = nnb
      mrs    = nnbf
      nrs    = nnbf
      nkl    = nnb
c
      ibix=0
      ibax=0
c
      lenax=nao(1)*nob(1)
      lenix=nco(1)*nob(1)
      if (lenax.gt.lbufso.or.lenix.gt.lbufso) then
         call lnkerr(' problems with lenax or lenix')
      end if
c
      mr1=mr-1
      ms1=ms-1
      mk1=mk-1
      ml1=ml-1
c
      if(itype.gt.2)go to 1000
      if(mr.ne.mmr.or.ms.ne.mms)go to 1000
c
c-----------------------------
c  locate the density matrices
c-----------------------------
c
c     write(iout,9010) mk,ml,mr,ms
c9010 format(' mcdnj mk ml mr ms ',4i5)
c
      mjk=0
      natest=nao(mk)*nao(ml)*nao(mr)*nao(ms)
      if(natest.eq.0)go to 411
c
      call mclden(ipqrs,mk,ml,nd,nstep,nij,mjk,
     $     ntot,ittype,nu,isqr,nblkd)
c
c
      if(idtype.eq.ittype)go to 411
      write(iout,409)
 409  format('  mcdnj  idtype ne ittype .. stop ')
      call lnkerr(' ')
 411  continue
c
c
c------------------------------------------c
c     idtype  labels an integral block
c             = 1    1 1 1 1
c             = 2    1 2 1 2
c             = 3    1 1 2 2
c             = 4    1 2 3 4
c                    1 1 2 3
c                    1 2 3 3
c-------------------------------------------c
c
      icrap=0
      ikrap=0
c
c
      ixab=0
      iabc=0
      if(noci.ne.0)   go to 100
      if(mk.ne.ml)    go to 50
      if(mr.ne.ms)    go to 50
      if(nao(mk).eq.0)go to 50
      if(nco(mr).eq.0)go to 50
      ncv=nao(mr)+nvirt(mr)
      if(ncv.eq.0)    go to 50
      ixab=1
      factiv=four*deg(mr)
      n34=(nao(mk)*(nao(mk)+1))/2
      n12=(nco(mr)*(nco(mr)+1))/2
c
cc    write(nfiv) mll,mk1,ml1,mr1,ms1,nco(mr),nob(mr),nco(ms),nob(ms),
cc   1 nao(mk),nao(ml),n12,n34,nblkd
c
      ihd(1)=mll
      ihd(2)=mk1
      ihd(3)=ml1
      ihd(4)=mr1
      ihd(5)=ms1
      ihd(6)=nco(mr)
      ihd(7)=nob(mr)
      ihd(8)=nco(ms)
      ihd(9)=nob(ms)
      ihd(10)=nao(mk)
      ihd(11)=nao(ml)
      ihd(12)=n12
      ihd(13)=n34
      ihd(14)=nblkd
c
c  abix will reside on rwf after aibx is added to abix in routine
c  mcdnk .. abix is temporarily stored on mcscr
c
      call iosys('write integer "header abix" to rwf',14,ihd,0,0)
c
      call iosys('does abix exist on mcscr',0,0,0,ians)
      if (ians.eq.'no') then
         ibnr=nco(1)*nob(1)
         ibnnk=(nao(1)+1)*nao(1)/2
         lpass=lbufso/ibnr
         npass=(ibnnk-1)/lpass+1
c
         call iosys('create real abix on mcscr',npass*lbufso,0,0,0)
      end if
      call iosys('rewind abix on mcscr',0,0,0,0)
c51
c
 50   continue
      if(natest.eq.0) go to 100
      iabc=1
      if(mr.eq.ms) go to 60
      n12=nao(mr)*nao(ms)
      go to 70
 60   continue
      n12=(nao(mr)*(nao(mr)+1))/2
 70   continue
      if(mk.eq.ml) go to 61
      n34=nao(mk)*nao(ml)
      go to 71
 61   continue
      n34=(nao(mk)*(nao(mk)+1))/2
 71   continue
cc    write(nfav) mll,mr1,ms1,mk1,ml1,nao(mr),nob(mr),nao(ms),nob(ms),
cc   1 nao(mk),nao(ml),n12,n34,nblkd
      ihd(1)=mll
      ihd(2)=mr1
      ihd(3)=ms1
      ihd(4)=mk1
      ihd(5)=ml1
      ihd(6)=nao(mr)
      ihd(7)=nob(mr)
      ihd(8)=nao(ms)
      ihd(9)=nob(ms)
      ihd(10)=nao(mk)
      ihd(11)=nao(ml)
      ihd(12)=n12
      ihd(13)=n34
      ihd(14)=nblkd
c52      call swrit(nfav,ihd,14)
      call iosys('write integer "header abcx" to rwf',14,ihd,0,0)
      call iosys('does abcx exist on rwf',0,0,0,ians)
      if (ians.eq.'no') then
         ibnr=nao(1)*nob(1)
         ibnnk=(nao(1)+1)*nao(1)/2
         lpass=lbufso/ibnr
         npass=(ibnnk-1)/lpass+1
c
         call iosys('create real abcx on rwf',npass*lbufso,0,0,0)
      else
         call iosys('rewind abcx on rwf',0,0,0,0)
      end if
 100  continue
c
      naol=nao(ml)
      naok=nao(mk)
      ncol=nco(ml)
      ncok=nco(mk)
      naor=nao(mr)
      naos=nao(ms)
      ncor=nco(mr)
      ncos=nco(ms)
      nvirtr=nvirt(mr)
      nvirts=nvirt(ms)
      nobr=nob(mr)
      nobs=nob(ms)
      ipr=isymc(mr)
      ips=isymc(ms)
      iread=0
      if(mr.ne.ms) go to 220
      lhaa=1
      lhac=((naor*(naor+1))/2          )*mrs+1
      lhcc=((naor*(naor+1))/2+naor*ncor)*mrs+1
      go to 230
 220  continue
      lhaa=1
      lhac=lhaa+(naor*naos)*mrs
      lhca=lhac+(naor*ncos)*mrs
      lhcc=lhca+(ncor*naos)*mrs
 230  continue
c
      if(mk.ne.ml)go to 600
c--------------------------
c   diagonal symmetry block
c--------------------------
      id=0
      nstart=0
c----------------------------------------c
c    read in mkl  coulomb matrices       c
c----------------------------------------c
      if(jkcore.eq.0)then
         iijt=nbufj/nnp
      else
         iijt=jmkl
         ncolj=(jnij-1)/jmij+1
         nbufj=ncolj*jmij*jmkl
      endif
      jst=jstrt
c
      lxj=1
      call mcrdjk(xjbuf,nbufj,ndab,jst,jkcore)
c
      do 500 k=1,nok
         do 495 l=1,k
            iread=iread+1
            if(iread.le.iijt) go to 405
            lxj=1
            call mcrdjk(xjbuf,nbufj,ndab,jst,jkcore)
            iread=1
 405        continue
            if(jkcore.eq.0) then
               call trtosq(xj,xjbuf(lxj),nbf,nnp)
            else
               call scopy(jnij,xjbuf(iread),jmkl,temp,1)
               call trtosq(xj,temp,nbf,nnp)
            endif
c
            ka=noba(lnoba(mk)+k)
            la=noba(lnoba(ml)+l)
            kala=ka+la
            if(kala.eq.0)go to 480
            if(ka.eq.0.or.la.eq.0)go to 470
            if(natest.eq.0)go to 425
            if(iabc.eq.0) go to 416
            if(mr.ne.ms) go to 410
            isym=0
            call mcabcx(c(ips),c(ipr),xj,nbr,nbr,naor,naor,
     $           nobr,nobr,temp,tv,nfav,isym,ncor,ncor,bufax(ibax+1))
            go to 415
 410        continue
            isym=1
            call mcabcx(c(ipr),c(ips),xj,nbr,nbs,naor,naos,
     $           nobr,nobs,temp,tv,nfav,isym,ncor,ncos,bufax(ibax+1))
 415        continue
            ibax=ibax+lenax
            if (ibax+lenax.gt.lbufso) then
               call iosys('write real abcx to rwf without rewinding',
     $              lbufso,bufax,0,0)
               ibax=0
            end if
 416        continue
c
c
c
c
            if(nstep.ne.1)go to 420
c
            call mcdxaa(den(id+nd),xj,nij,mrs,nstep,itran,thrsh,
     $           hess(lhaa),nstart,ntot)
            id=id+nij
c
            go to 425
c
 420        continue
c
            nstart=nstart+1
ccc
            call mcdxaa(den(nd),xj,nij,mrs,nstep,itran,thrsh,
     $           hess(lhaa),nstart,ntot)
ccc
c
 425        continue
            if(mr.ne.ms)go to 490
            if(nu.gt.1) go to 490
            if(ncor.eq.0)go to 490
c-------------------------------------------
c   xj contributes to the fab  fock-operator
c-------------------------------------------
            ldab1=ldab(mk)+(ka-1)*naol+la-1
            factj=two*deg(mr)*dab(ldab1)
            if(k.eq.l) factj=factj*pt5
            lfab1=lfab(mr)
c
            call mcstv(fab(lfab1),xj,factj,mrs)
c
            go to 490
 470        continue
c
c
            if(idtype.gt.1) go to 490
c
c---------------------
c    core-active terms
c---------------------
c
            lden1=ldab(mk)+(ka-1)*naok
            lh=lhac+(l-1)*naok*mrs
c
            call mcjac(dab(lden1),hess(lh),xj,naok,mrs,itran,thrsh)
c
            go to 490
c
 480        continue
c-------------------
c    core-core terms
c-------------------
c
            if(idtype.ne.1) go to 490
            if(k.ne.l)go to 485
c
            lh=lhcc+((k*(k+1))/2)*mrs-mrs
            call mcstv(hess(lh),xj,aii(mk),mrs)
c
            go to 490
 485        continue
c
            lh=lhcc+((k*(k-1))/2)*mrs+(l-1)*mrs
            call mcstv(hess(lh),xj,aij(mk,mk),mrs)
c
 490        continue
c--------------------------------------
c   core contribution to abii  i-vector
c--------------------------------------
            if(ixab.eq.0)go to 495
            if(ka.eq.0.or.la.eq.0)go to 495
            call mcabix(c(ipr),xj,nbr,ncor,nobr,temp,tv,nfiv,
     $           factiv,ikt,'abix',bufix(ibix+1))
            ibix=ibix+lenix
            if (ibix+lenix.gt.lbufso) then
               call iosys('write real abix to mcscr without rewinding',
     $              lbufso,bufix,0,0)
               ibix=0
            end if
c
c
 495     lxj=lxj+nnp
 500  continue
      go to 1000
 600  continue
c------------------------------
c   off-diagonal symmetry block
c------------------------------
      id=0
      nstart=0
c
c----------------------------------------c
c    read in a block coulomb integrals   c
c----------------------------------------c
c
      iijt=nbufj/nrs
      lxj=1
      call mcrdjk(xjbuf,nbufj,ndab)
      do 650 l=1,nol
         la=noba(lnoba(ml)+l)
         do 640 k=1,nok
            ka=noba(lnoba(mk)+k)
            iread=iread+1
            if(iread.lt.iijt) go to 601
            lxj=1
            call mcrdjk(xjbuf,nbufj,ndab)
 601        continue
c
c
            kala=ka+la
            if(kala.eq.0)go to 630
            if(ka.eq.0.or.la.eq.0)go to 610
c---------------------
c  active-active block
c---------------------
c
            if(natest.eq.0)go to 640
            if(iabc.eq.0) go to 716
            if(mr.ne.ms) go to 710
            isym=0
            call mcabcx(c(ips),c(ipr),xj,nbr,nbr,naor,naor,
     $           nobr,nobr,temp,tv,nfav,isym,ncor,ncor,bufax(ibax+1))
            go to 715
 710        continue
            isym=1
            call mcabcx(c(ipr),c(ips),xj,nbr,nbs,naor,naos,
     $           nobr,nobs,temp,tv,nfav,isym,ncor,ncos,bufax(ibax+1))
 715        continue
            ibax=ibax+lenax
            if (ibax+lenax.gt.lbufso) then
               call iosys('write real abcx to rwf without rewinding',
     $              lbufso,bufax,0,0)
               ibax=0
            end if
 716        continue
c
c
            if(nstep.ne.1)go to 605
c
            call mcdxaa(den(id+nd),xj,nij,mrs,nstep,itran,thrsh,
     $           hess(lhaa),nstart,ntot)
c
            id=id+nij
            go to 640
 605        continue
c
            nstart=nstart+1
            call mcdxaa(den(nd),xj,nij,mrs,nstep,itran,thrsh,
     $           hess(lhaa),nstart,ntot)
c
            go to 640
 610        continue
            if(idtype.gt.2)go to 640
c-------------------
c  active-core block
c-------------------
            if(ka.eq.0)go to 620
            lden=ldab(mk)+(ka-1)*naok
            lh=lhac+(l-1)*naok*mrs
            call mcjac (dab(lden),hess(lh),xj,naok,mrs,itran,thrsh)
c
            go to 640
c
 620        continue
            lden=ldab(ml)+(la-1)*naol
            lh=lhca+(k-1)*naol*mrs
            call mcjac (dab(lden),hess(lh),xj,naol,mrs,itran,thrsh)
            go to 640
c
 630        continue
c-------------------
c  core-core mixings
c-------------------
            if(idtype.gt.2)go to 640
            lh=lhcc+(l-1)*ncok*mrs+(k-1)*mrs
            call mcstv(hess(lh),xj,aij(mk,ml),mrs)
c
 640     lxj=lxj+nrs
 650  continue
 1000 continue
c
      if (ibax.ne.0) then
         call iosys('write real abcx to rwf without rewinding',
     $        lbufso,bufax,0,0)
      end if
      if (ibix.ne.0) then
         call iosys('write real abix to mcscr without rewinding',
     $        lbufso,bufix,0,0)
      end if
c
      return
      end
