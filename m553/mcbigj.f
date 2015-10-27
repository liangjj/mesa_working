*deck @(#)mcbigj.f	5.1  11/6/94
      subroutine mcbigj(xj,hess,den,aii,aij,
     $     dab,ldab,deg,fab,lfab,noba,lnoba,nco,nao,nvirt,nob,c,isymc,
     $     nsym,noci,temp,tv,itran,thrsh,mmr,mms,ndab,nfiv,nfav,
     $     lhaa,lhac,lhca,lhcc,
     $     xjbuf,nocc,nbf,nbufj,jstrt,jnij,jmij,jnkl,jmkl,jkcore,
     $     bufix,bufax,lbufso,lxj,iread,jstsav,ird,nbufjj,mblkd,lblkd)
c
c***begin prologue     mcbigj
c***date written       871022   (yymmdd)
c***revision date      900417   (yymmdd)
c
c   17 april,   1990   rlm at lanl
c     passing mblkd,lblkd to avoid mclden as an entry point.
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcbigj.f	5.1   11/6/94
c
c***purpose
c
c***description
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
c***references
c
c***routines called    (none)
c
c***end prologue       mcbigj
c
      implicit real*8 (a-h,o-z)
c
      dimension xj(2),hess(2),den(2),dab(2),ldab(2),deg(2),fab(2),
     $     lfab(2),nco(2),nao(2),nvirt(2),nob(2),c(2),isymc(2),temp(2),
     $     tv(2),noba(2),lnoba(2),aij(10,2),aii(2),xjbuf(2)
      integer mblkd(51,2),lblkd(2)
      real*8 bufix(lbufso),bufax(lbufso)
c
      character*3 ians
      dimension ihd(14)
c
      common /io/ inp,iout
      common / number / zero,pt5,one,two,four,eight
c
      ikt=0
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
         call lnkerr('intrn:  problems with lenax or lenix')
      end if
c
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
      mjk=0
c
      call mclden(ipqrs,mk,ml,nd,nstep,nij,mjk,
     $     ntot,ittype,nu,isqr,nblkd,mblkd,lblkd)
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
      call iosys('write integer "header abix" to rwf',14,ihd,0,' ')
c
      call iosys('does abix exist on mcscr',0,0,0,ians)
      if (ians.eq.'no') then
         ibnr=nco(1)*nob(1)
         ibnnk=(nao(1)+1)*nao(1)/2
         lpass=lbufso/ibnr
         npass=(ibnnk-1)/lpass+1
c
         call iosys('create real abix on mcscr',npass*lbufso,0,0,' ')
      end if
      call iosys('rewind abix on mcscr',0,0,0,' ')
c51
c
 50   continue
c
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
c
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
      call iosys('write integer "header abcx" to rwf',14,ihd,0,' ')
      call iosys('does abcx exist on rwf',0,0,0,ians)
c
      if (ians.eq.'no') then
         ibnr=nao(1)*nob(1)
         ibnnk=(nao(1)+1)*nao(1)/2
         lpass=lbufso/ibnr
         npass=(ibnnk-1)/lpass+1
c
         call iosys('create real abcx on rwf',npass*lbufso,0,0,' ')
      else
         call iosys('rewind abcx on rwf',0,0,0,' ')
      end if
c
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
      lhaa=1
c
c--------------------------
c   diagonal symmetry block
c--------------------------
      id=0
      nstart=0
c
c----------------------------------------c
c    read in mkl  coulomb matrices       c
c----------------------------------------c
c        nbf2=nbf*(nbf+1)/2
         noc2=nocc*(nocc+1)/2
         iijt=nbufj/nnp
         nbufj=iijt*nnp
         nbufjj=nbufj
         if(iijt.gt.noc2) then
          iijt=noc2
          nbufjj=noc2*nnp
          nbufj=nbufjj
         endif
         ird=iijt
c
c
      lxj=1
      jst=jstrt
      jstsav=jst
c
      call mcrdj(xjbuf,nbufjj,ndab,jst,jkcore)
c
c
      do 500 k=1,naok
         do 495 l=1,k
            iread=iread+1
            if(iread.le.iijt) go to 405
            ileft=min(iijt,noc2-ird)
            ird=ird+ileft
            nbufjj=ileft*nnp
            lxj=1
            jstsav=jst
            call mcrdj(xjbuf,nbufjj,ndab,jst,jkcore)
            iread=1
 405        continue
c
               call trtosq(xj,xjbuf(lxj),nbf,nnp)
c
c
            if(iabc.eq.0) go to 416
            isym=0
            call mcabcx(c(ips),c(ipr),xj,nbr,nbr,naor,naor,
     $           nobr,nobr,temp,tv,nfav,isym,ncor,ncor,bufax(ibax+1))
            ibax=ibax+lenax
            if (ibax+lenax.gt.lbufso) then
               call iosys('write real abcx to rwf without rewinding',
     $              lbufso,bufax,0,' ')
               ibax=0
            end if
 416        continue
c
c
c
c-2      write(iout,*)'  nstep  ',nstep
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
c
c-------------------------------------------
c   xj contributes to the fab  fock-operator
c-------------------------------------------
c
            ldab1=ldab(mk)+(k-1)*naol+l-1
            factj=two*deg(mr)*dab(ldab1)
            if(k.eq.l) factj=factj*pt5
            lfab1=lfab(mr)
c
            call mcstv(fab(lfab1),xj,factj,mrs)
c
 490        continue
c--------------------------------------
c   core contribution to abii  i-vector
c--------------------------------------
            if(ixab.eq.0)go to 495
            call mcabix(c(ipr),xj,nbr,ncor,nobr,temp,tv,nfiv,
     $           factiv,ikt,'abix',bufix(ibix+1))
            ibix=ibix+lenix
            if (ibix+lenix.gt.lbufso) then
               call iosys('write real abix to mcscr without rewinding',
     $              lbufso,bufix,0,' ')
               ibix=0
            end if
c
c
 495     lxj=lxj+nnp
 500  continue
 1000 continue
c
      if (ibax.ne.0) then
         call iosys('write real abcx to rwf without rewinding',
     $        lbufso,bufax,0,' ')
      end if
      if (ibix.ne.0) then
         call iosys('write real abix to mcscr without rewinding',
     $        lbufso,bufix,0,' ')
      end if
c
      return
      end
