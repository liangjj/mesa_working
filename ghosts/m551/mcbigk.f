*deck  %W% %G%
      subroutine mcbigk(txk,xk,hess,den,bii,bij,
     $     dab,ldab,deg,fab,lfab,noba,lnoba,nco,nao,nvirt,nob,c,isymc,
     $     nsym,noci,temp,tv,itran,thrsh,xkt,bijt,mmr,mms,ndab,
     $     nfiv,lhaa,lhac,lhca,lhcc,
     $     nocc,nbf,nbufk, kstrt,knij,kmij,knkl,kmkl,jkcore,
     $     bufix,bufab,lbufso,lxk,iread,kstsav,ird,nbufkk,mblkd,lblkd)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      900417   (yymmdd)
c
c   april 17,   1990   rlm at lanl
c     passing mblkd,lblkd to avoid mclden as an entry point.
c***keywords
c***author             lengsfield, byron (brl)
c***source              %W% %G%
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
cc
      implicit real*8 (a-h,o-z)
      dimension xk(2),hess(2),den(2),dab(2),ldab(2),deg(2),fab(2),
     $     lfab(2),nco(2),nao(2),nvirt(2),c(2),isymc(2),temp(2),tv(2),
     $     noba(2),lnoba(2),bij(10,2),bii(2),
     $     xkt(2),bijt(10,2),nob(2),txk(2)
      integer mblkd(51,2), lblkd(2)
      real*8 bufix(lbufso),bufab(lbufso)
c
      character*3 ians
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
c
c-----------------------------
c  locate the density matrices
c-----------------------------
c
c
      mjk=1
c
      call mclden(ipqrs,mk,ml,nd,nstep,nij,mjk,
     $     ntot,ittype,nu,isqr,nblkd,mblkd,lblkd)
c
      if(ittype.eq.idtype)go to 411
      write(iout,409)
 409  format('  mcdnk  idtype ne ittype .. stop')
      call lnkerr(' ')
 411  continue
c
c
      nleft=0
c
      ipr=isymc(mr)
      ips=isymc(ms)
      ixab=0
      if(noci.ne.0)   go to 50
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
         call iosys('create real abix on rwf',npass*lbufso,0,0,' ')
      end if
      call iosys('rewind abix on mcscr',0,0,0,' ')
      call iosys('rewind abix on rwf',0,0,0,' ')
      call iosys('read real abix from mcscr without rewinding',
     $     lbufso,bufab,0,' ')
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
      lhaa=1
c
c--------------------------
c   diagonal symmetry block
c--------------------------
c
      id=0
      idt=0
c
      nstart=0
      nstrtt=0
c
         noc2=nocc*(nocc+1)/2
         iijt=nbufk/nrs
         nbufk=iijt*nrs
         nbufkk=nbufk
c
         if(iijt.gt.noc2) then
          iijt=noc2
          nbufkk=noc2*nrs
          nbufk=nbufkk
         endif
c
      lxk=1
      kst=kstrt
      kstsav=kst
      call mcrdk(xk,nbufkk,ndab,kst,jkcore)
      ird=iijt
c
      do 500 k=1,naok
         do 495 l=1,k
            iread=iread+1
            if(iread.le.iijt) go to 405
            ileft=min(iijt,noc2-ird)
            nbufkk=ileft*nrs
            ird=ird+ileft
            lxk=1
            iread=1
            kstsav=kst
            call mcrdk(xk,nbufkk,ndab,kst,jkcore)
 405        continue
c
c
            if(nstep.ne.1)go to 415
c
c------------------------------------------------------
c   multiply the exchange operator times the appropiate
c   density matrix elements
c------------------------------------------------------
c
            id =((l-1)*naol+k-1)*nij
            call mcdxaa(den(id+nd),xk(lxk),nij,mrs,nstep,itran,thrsh,
     $           hess(lhaa),nstart,ntot)
c
            if(k.eq.l)   go to 410
c
c---------------------------------------------------------
c   the transpose of this exchange operator is also needed
c---------------------------------------------------------
c
            idt=((k-1)*naol+l-1)*nij
            call mctrsp(xk(lxk),xkt,nbr,nbs)
            call mcdxaa(den(idt+nd),xkt,nij,mrs,nstep,itran,thrsh,
     $           hess(lhaa),nstrtt,ntott)
cc
 410        continue
cc
cc
            go to 425
 415        continue
            write(iout,416)
 416        format(' problem mcdnk transpose d.m. ? ')
            call lnkerr(' ')
c
 425        continue
c
c-------------------------------------------
c   xk contributes to the fab  fock-operator
c-------------------------------------------
c
            ldab1=ldab(mk)+(k-1)*naol+l-1
            lfab1=lfab(mr)
            if(k.eq.l)go to 455
c
            factk=-deg(mr)*dab(ldab1)
c
            call mcstv(fab(lfab1),xk(lxk),factk,mrs)
            go to 490
c
 455        continue
c
            factk=-pt5*deg(mr)*dab(ldab1)
            call mcstv(fab(lfab1), xk(lxk),factk,mrs)
c
 490        continue
c
c--------------------------------------
c   core contribution to abii  i-vector
c--------------------------------------
c
            if(ixab.eq.0)go to 495
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
     $              lbufso,bufab,0,' ')
               if(kpass.lt.npass) then
                  call iosys('read real abix from mcscr '//
     $                 'without rewinding',lbufso,bufab,0,' ')
                  kpass=kpass+1
               end if
               ibix=0
            end if
c
c
 495     lxk=lxk+nrs
 500  continue
c
c
      if (ibix.ne.0) then
         call iosys('write real abix to rwf without rewinding',
     $        lbufso,bufab,0,' ')
      end if
      return
      end
