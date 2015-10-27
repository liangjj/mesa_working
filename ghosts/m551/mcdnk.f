*deck %W%  %G%
      subroutine mcdnk(txk,xk,hess,den,bii,bij,
     $     dab,ldab,deg,fab,lfab,noba,lnoba,nco,nao,nvirt,nob,c,isymc,
     $     nsym,noci,temp,tv,itran,thrsh,xkt,bijt,mmr,mms,ndab,
     $     nfiv,lhaa,lhac,lhca,lhcc,
     $     nocc,nbf,nbufk, kstrt,knij,kmij,knkl,kmkl,jkcore,
     $     bufix,bufab,lbufso)
C
C***Begin prologue
C***Date written       871022   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C
C***Keywords
C***Author             Lengsfield, Byron (BRL)
C***Source             %W%   %G%
C
C***Purpose
C
C***Description
C
C***References
C
C***Routines called    (none)
C
C***End prologue
C
      implicit real*8(a-h,o-z)
cc
      dimension xk(2),hess(2),den(2),dab(2),ldab(2),deg(2),fab(2),
     $     lfab(2),nco(2),nao(2),nvirt(2),c(2),isymc(2),temp(2),tv(2),
     $     noba(2),lnoba(2),bij(10,2),bii(2),mklrs(20,2),mrskl(2),
     $     xkt(2),bijt(10,2),nob(2),txk(2)
      real*8 bufix(lbufso),bufab(lbufso)
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
c  this program processes the exchange integral blocks
c  these blocks contribute to the orbital hessian and
c  to the  fab  fock operator. in addition, some of these
c  blocks are needed to compute derivatives of fock matrix
c  elements. these last terms are needed when ci-coupling
c  is invoked and are called i-vector terms in our convention
c
c------------------------------------------------------------
c      do 1000 irskl=1,nklv
c
c     irskl=mrskl(iklv)
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
      itype  = 3
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
      lenix=nco(1)*nob(1)
      if (lenix.gt.lbufso) call lnkerr('problem with lenix')
      mr1=mr-1
      ms1=ms-1
      mk1=mk-1
      ml1=ml-1
c
      if(itype.lt.3)go to 1000
      if(mr.ne.mmr.or.ms.ne.mms)go to 1000
c-----------------------------
c  locate the density matrices
c-----------------------------
c
c     write(iout,9010) mk,ml,mr,ms
c9010 format('  mcdnk  mk ml mr ms ',4i5)
c
      mjk=1
      natest=nao(mr)*nao(ms)*nao(mk)*nao(ml)
      if(natest.eq.0)go to 411
c
      call mclden(ipqrs,mk,ml,nd,nstep,nij,mjk,
     $     ntot,ittype,nu,isqr,nblkd)
c
      if(ittype.eq.idtype)go to 411
      write(iout,409)
 409  format('  mcdnk  idtype ne ittype .. stop')
      call lnkerr(' ')
 411  continue
c
c
      idsk1=idad0
      nleft=0
c
      ipr=isymc(mr)
      ips=isymc(ms)
      ixab=0
      if(noci.ne.0)   go to 50
      if(mk.ne.ml)    go to 50
      if(mr.ne.ms)    go to 50
      if(nao(mk).eq.0)go to 50
      if(nco(mr).eq.0)go to 50
      ncv=nao(mr)+nvirt(mr)
      if(ncv.eq.0)    go to 50
      ixab=1
      fctiv1=-two*deg(mr)
      fctiv2=-one*deg(mr)
      n34=(nao(mk)*(nao(mk)+1))/2
      n12=(nco(mr)*(nco(mr)+1))/2
      call iosys('does abix exist on rwf',0,0,0,ians)
      ibnr=nco(1)*nob(1)
      ibnnk=(nao(1)+1)*nao(1)/2
      lpass=lbufso/ibnr
      npass=(ibnnk-1)/lpass+1
c
      if (ians.eq.'no') then
         call iosys('create real abix on rwf',npass*lbufso,0,0,0)
      end if
      call iosys('rewind abix on mcscr',0,0,0,0)
      call iosys('rewind abix on rwf',0,0,0,0)
      call iosys('read real abix from mcscr without rewinding',
     $     lbufso,bufab,0,0)
      kpass=1
c
 50   continue
c
      naol=nao(ml)
      naok=nao(mk)
      ncol=nco(ml)
      ncok=nco(mk)
      naor=nao(mr)
      naos=nao(ms)
      ncor=nco(mr)
      ncos=nco(ms)
      iread=0
c
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
      if(mk.ne.ml)go to 600
c--------------------------
c   diagonal symmetry block
c--------------------------
c
c
      id=0
      idt=0
c
      nstart=0
      nstrtt=0
      if(jkcore.eq.0) then
         iijt=nbufk/nrs
      else
         write(iout,*)' bug in mcdnk  jkcore ne 0  stop '
         call lnkerr(' ')
      endif
c
      lxk=1
      kst=kstrt
      call mcrdjk(xk,nbufk,ndab,kst,jkcore)
c
      do 500 k=1,nok
         do 495 l=1,k
            iread=iread+1
            if(iread.le.iijt) go to 405
            lxk=1
            iread=1
            call mcrdjk(xk,nbufk,ndab,kst,jkcore)
 405        continue
c
            lhxx=lhaa+mrs
c
            ka=noba(lnoba(mk)+k)
            la=noba(lnoba(ml)+l)
            kala=ka+la
            if(kala.eq.0)go to 480
            if(ka.eq.0.or.la.eq.0)go to 470
            if(natest.eq.0)go to 425
c
            if(nstep.ne.1)go to 415
c
c------------------------------------------------------
c   multiply the exchange operator times the appropiate
c   density matrix elements
c------------------------------------------------------
c
            id =((la-1)*naol+ka-1)*nij
            call mcdxaa(den(id+nd),xk(lxk),nij,mrs,nstep,itran,thrsh,
     $           hess(lhaa),nstart,ntot)
c
            if(k.eq.l)   go to 410
c
c---------------------------------------------------------
c   the transpose of this exchange operator is also needed
c---------------------------------------------------------
c
            idt=((ka-1)*naol+la-1)*nij
            call mctrsp(xk(lxk),xkt,nbr,nbs)
            call mcdxaa(den(idt+nd),xkt,nij,mrs,nstep,itran,thrsh,
     $           hess(lhaa),nstrtt,ntott)
cc
 410        continue
c
            go to 425
 415        continue
            write(iout,416)
 416        format(' problem mcdnk transpose d.m. ? ')
            call lnkerr(' ')
c
 425        continue
            if(mr.ne.ms)go to 490
c
c-------------------------------------------
c   xk contributes to the fab  fock-operator
c-------------------------------------------
c
            ldab1=ldab(mk)+(ka-1)*naol+la-1
            lfab1=lfab(mr)
            if(k.eq.l)go to 455
c
            factk=-deg(mr)*dab(ldab1)
            call mcstv(fab(lfab1),xk(lxk),factk,mrs)
            go to 490
c
 455        continue
c
            factk=-pt5*deg(mr)*dab(ldab1)
            call mcstv(fab(lfab1), xk(lxk),factk,mrs)
c
            go to 490
 470        continue
c
c---------------------
c    core-active terms
c---------------------
c
            if(idtype.gt.1)go to 490
            lden1=ldab(mk)+(ka-1)*naok
            lh=lhac+(l-1)*naok*mrs
            factk=four*deg(mk)
c
            call mckac(dab(lden1),hess(lh),xk(lxk),naok,mrs,itran,
     $           thrsh,factk)
            factk=-deg(mk)
            call mcktac(dab(lden1),hess(lh),xk(lxk),naok,nbr,nbs,itran,
     $           thrsh,factk)
c
            go to 490
c
 480        continue
c
c-------------------
c    core-core terms
c-------------------
c
            if(idtype.ne.1) go to 490
            if(k.ne.l)go to 485
c
            lh=lhcc+((k*(k+1))/2)*mrs-mrs
            call mcstv(hess(lh),xk(lxk),bii(mk),mrs)
c
            go to 490
c
 485        continue
c
            lh=lhcc+((k*(k-1))/2)*mrs+(l-1)*mrs
            call mcstv(hess(lh), xk(lxk), bij(mk,mk),mrs)
            call mcstvt(hess(lh),xk(lxk),bijt(mk,mk),nbr,nbs)
c
 490        continue
c
c--------------------------------------
c   core contribution to abii  i-vector
c--------------------------------------
c
            if(ixab.eq.0)go to 495
            if(ka.eq.0.or.la.eq.0)go to 495
            if(k.eq.l) go to 491
            ikt=1
            call mcaibx(c(ipr),xk(lxk),nbr,nco(mr),nob(mr),temp,tv,nfiv,
     $           fctiv2,ikt,'aibx',bufix(ibix+1),bufab(ibix+1))
            go to 494
 491        continue
            ikt=0
            call mcaibx(c(ipr),xk(lxk),nbr,nco(mr),nob(mr),temp,tv,
     $           nfiv,fctiv1,ikt,'aibx',bufix(ibix+1),bufab(ibix+1))
c
 494        continue
            ibix=ibix+lenix
            if (ibix+lenix.gt.lbufso) then
               call iosys('write real abix to rwf without rewinding',
     $              lbufso,bufab,0,0)
               if(kpass.lt.npass) then
                  call iosys('read real abix from mcscr '//
     $                 'without rewinding',lbufso,bufab,0,0)
                  kpass=kpass+1
               end if
               ibix=0
            end if
c
c
 495     lxk=lxk+nrs
 500  continue
      go to 1000
 600  continue
c
c------------------------------
c   off-diagonal symmetry block
c------------------------------
c
      id=0
      idt=0
      nstart=0
c
      iijt=nbufk/nrs
      lxk=1
      call mcrdjk(xk,nbufk,ndab,kst,jkcore)
c
      do 650 l=1,nol
         la=noba(lnoba(ml)+l)
         do 640 k=1,nok
            ka=noba(lnoba(mk)+k)
            iread=iread+1
            if(iread.lt.iijt) go to 601
            lxk=1
            call mcrdjk(xk,nbufk,ndab,kst,jkcore)
 601        continue
c
c
            kala=ka+la
            if(kala.eq.0)go to 630
            if(ka.eq.0.or.la.eq.0)go to 610
c
c---------------------
c  active-active block
c---------------------
c
            if(natest.eq.0)go to 640
            if(nstep.ne.1)go to 605
cc
            call mcdxaa(den(id+nd),xk(lxk),nij,mrs,nstep,itran,thrsh,
     $           hess(lhaa),nstart,ntot)
            id=id+nij
cc
            go to 640
 605        continue
c
c
cc
            nstart=nstart+1
            call mcdxaa(den(nd),xk(lxk),nij,mrs,nstep,itran,thrsh,
     $           hess(lhaa),nstart,ntot)
cc
            go to 640
 610        continue
            if(idtype.gt.3)go to 640
c-------------------
c  active-core block
c-------------------
c
            if(ka.eq.0)go to 620
            lden=ldab(mk)+(ka-1)*naok
            lh=lhac+(l-1)*naok*mrs
            if(idtype.eq.3)go to 615
            factk=-deg(ml)
            call mckac(dab(lden),hess(lh),xk(lxk),naok,mrs,itran,
     $           thrsh,factk)
            go to 640
 615        continue
            factk=four*deg(ml)
            call mckac(dab(lden),hess(lh),xk(lxk),naok,mrs,itran,thrsh,
     $           factk)
            go to 640
c
 620        continue
            lden=ldab(ml)+(la-1)*naol
            lh=lhca+(k-1)*naol*mrs
            if(idtype.eq.3)go to 625
            factk=-deg(mk)
            call mckac(dab(lden),hess(lh),xk(lxk),naol,mrs,itran,
     $           thrsh,factk)
            go to 640
 625        continue
            factk=four*deg(mk)
            call mckac(dab(lden),hess(lh),xk(lxk),naol,mrs,itran,thrsh,
     $           factk)
            go to 640
c
 630        continue
c-------------------
c  core-core mixings
c-------------------
            if(idtype.gt.3)go to 640
            lh=lhcc+((l-1)*ncok+k)*mrs-mrs
            if(idtype.eq.3)go to 635
c-------------------------c
c   ijij symmetry block
c-------------------------c
            call mcstv(hess(lh),xk(lxk),bijt(mk,ml),mrs)
            go to 640
 635        continue
c-------------------------c
c   iijj symmetry block
c-------------------------c
            call mcstv(hess(lh),xk(lxk), bij(mk,ml),mrs)
c
 640     lxk=lxk+nrs
 650  continue
 1000 continue
c
c
      if (ibix.ne.0) then
         call iosys('write real abix to rwf without rewinding',
     $        lbufso,bufab,0,0)
      end if
      return
      end
