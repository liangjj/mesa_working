*deck @(#)mccor.f	5.1  11/6/94
      subroutine mccor(xj,txk,xk,hess,den,aii,aij,bii,bij,bijt,
     $     dab,ldab,deg,fab,lfab,noba,lnoba,nco,nao,nvirt,nob,c,isymc,
     $     nsym,noci,temp,tv,itran,thrsh,mmr,mms,ndab,nfiv,nfav,
     $     lhaa,lhac,lhca,lhcc,
     $     xjbuf,nocc,nbf,
     $     nbufj,jstrt,jnij,jmij,jnkl,jmkl,jkcore,
     $     nbufk,kstrt,knij,kmij,knkl,kmkl,
     $     lxjsav,lxksav,ireadj,ireadk,nbufcc,nbufac,
     $     f,incorh, irdj,irdk,nbufjj,nbufkk)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      871116   (yymmdd)
c   16 november 1987   bhl      lanl
c   modified to skip loops  if nao=0
c   ireadj ireadk lxjsav lxksav defined for closed shell scf run
c   changes denoted by c..bhl..close shell scf
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mccor.f	5.1   11/6/94
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
c
      implicit real*8(a-h,o-z)
c
      dimension xj(2),hess(2),den(2),dab(2),ldab(2),deg(2),fab(2),
     $     lfab(2),nco(2),nao(2),nvirt(2),nob(2),c(2),isymc(2),temp(2),
     $     tv(2),noba(2),lnoba(2),aij(10,2),aii(2),xjbuf(2),f(2),
     $     txk(2),xk(2),bij(10,2),bijt(10,2),bii(2)
c
      character*3 ians
c
      common / number / zero,pt5,one,two,four,eight
c
      common /io/ inp,iout
c
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
cc
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
c
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
c
c
      call iosys('does hess_cc exist on rwf',0,0,0,ians)
c
      if(ians.eq.'no') then
         ncot=nco(1)*(nco(1)+1)/2
         call iosys('create real hess_cc on rwf',ncot*mrs,0,0,' ')
      endif
c
      call iosys('rewind hess_cc on rwf',0,0,0,' ')
c
c..bhl..closed shell scf
c
      if(naok.ne.0) then
      call iosys('does hess_ac exist on rwf',0,0,0,ians)
c
      if(ians.eq.'no') then
         nconao=nco(1)*nao(1)
         call iosys('create real hess_ac on rwf',nconao*mrs,0,0,' ')
      endif
c
      call iosys('rewind hess_ac on rwf',0,0,0,' ')
      endif
c
c----------------------------------------c
c    read in mkl  coulomb matrices       c
c----------------------------------------c
c
         noc2=nocc*(nocc+1)/2
         iijt=nbufj/nnp
         nbufj=iijt*nnp
         if(iijt.gt.noc2) then
          iijt=noc2
          nbufj=iijt*nnp
         endif
c
c..bhl..closed shell scf
c
      if (naol.eq.0) then
         ireadj=0
         ireadk=0
         lxjsav=1
         lxksav=1
      end if
c
      jst=jstrt
      lxj=lxjsav
c
c..bhl..closed shell scf
c
      if(naol.eq.0) then
       nbufjj=nbufj
       irdj=iijt
      endif
c
      call mcrdj(xjbuf,nbufjj,ndab,jst,jkcore)
c----------------------------------------c
c    read in mkl  exchange matrices      c
c----------------------------------------c
         iikt=nbufk/nrs
         nbufk=iikt*nrs
         if(iikt.gt.noc2) then
          iikt=noc2
          nbufk=iikt*nrs
         endif
c
c..bhl..closed..shell
c
         if(naol.eq.0) then
           irdk=iikt
           nbufkk=nbufk
         endif
c
      kst=kstrt
c
      lxk=lxksav
      call mcrdk(xk,nbufkk,ndab,kst,jkcore)
c
c..bhl..closed shell scf
c
      if(naok.ne.0) then
c
      nwac=nbufac*mrs
      nwcc=nbufcc*mrs
      lhac=1
      lhcc=lhac+nwac
      call rzero(hess(lhac),nwac)
      call rzero(hess(lhcc),nwcc)
c
      else
c
      nwac=0
      nwcc=nbufcc*mrs
      lhac=1
      lhcc=lhac+nwac
      call rzero(hess(lhcc),nwcc)
c
      endif
c
      llhac=lhac
      llhcc=lhcc
c
      icc=0
      iac=0
c
c-2      write(iout,*)' do loop  ncok  naok ',ncok,naok
c
      do 500 k=1,ncok
c
c..bhl..closed shell scf
c
         if(naok.eq.0) go to 480
c
         do 495 l=1,naok
c
c  j-read
c
            ireadj=ireadj+1
            if(ireadj.le.iijt) go to 404
            ileftj=min(iijt,noc2-irdj)
            nbufjj=ileftj*nnp
            irdj=irdj+ileftj
            lxj=1
            call mcrdj(xjbuf,nbufjj,ndab,jst,jkcore)
            ireadj=1
 404        continue
c
c
c   k-read
c
            ireadk=ireadk+1
            if(ireadk.le.iikt) go to 405
            ileftk=min(iikt,noc2-irdk)
            nbufkk=ileftk*nrs
            irdk=irdk+ileftk
            lxk=1
            ireadk=1
            call mcrdk(xk,nbufkk,ndab,kst,jkcore)
 405        continue
c
c
               call trtosq(xj,xjbuf(lxj),nbf,nnp)
c
c---------------------
c   core-active terms
c---------------------
c
            lden1=ldab(mk)+(l-1)*naok
c
            call mcjac(dab(lden1),hess(llhac),xj,naok,mrs,itran,thrsh)
c
c
            factk=four*deg(mk)

c
               call mckac(dab(lden1),hess(llhac),xk(lxk),naok,mrs,
     $              itran,thrsh,factk)
               factk=-deg(mk)
               call mcktac(dab(lden1),hess(llhac),xk(lxk),naok,nbr,
     $              nbs,itran,thrsh,factk)
c
c
            lxk=lxk+nrs
 495     lxj=lxj+nnp
c
         llhac=llhac+naok*mrs
         iac=iac+naok
         if(iac.eq.nbufac.and.incorh.eq.0) then
c-2        write(iout,*)'  write to hess_ac ',nwac,hess(lhac)
            call iosys('write real hess_ac on rwf without rewinding',
     $           nwac,hess(lhac),0,' ')
            call rzero(hess(lhac),nwac)
            iac=0
            llhac=lhac
         endif
c
c
 480     continue
c-------------------
c    core-core terms
c-------------------
c
         do 496 l=1,k
c
c  j-read
c
            ireadj=ireadj+1
            if(ireadj.le.iijt) go to 4405
            lxj=1
            ileftj=min(iijt,noc2-irdj)
            nbufjj=ileftj*nnp
            irdj=irdj+ileftj
            call mcrdj(xjbuf,nbufjj,ndab,jst,jkcore)
            ireadj=1
 4405       continue
c
c   k-read
c
            ireadk=ireadk+1
            if(ireadk.le.iikt) go to 4406
            lxk=1
            ileftk=min(iikt,noc2-irdk)
            irdk=irdk+ileftk
            nbufkk=ileftk*nrs
            ireadk=1
            call mcrdk(xk,nbufkk,ndab,kst,jkcore)
 4406       continue
c
               call trtosq(xj,xjbuf(lxj),nbf,nnp)
c
            if(k.ne.l)go to 485
c
            call mcstv(hess(llhcc),xj,aii(mk),mrs)
c
               call mcstv(hess(llhcc),xk(lxk),bii(mk),mrs)
c
            call vadd(hess(llhcc),hess(llhcc),fab(lfab(mr)),mrs)
c
            call saxpy(mrs,deg(mr),f(lfab(mr)),1,hess(llhcc),1)
c
            go to 490
 485        continue
c
c
            call mcstv(hess(llhcc),xj,aij(mk,mk),mrs)
c
               call mcstv(hess(llhcc), xk(lxk), bij(mk,mk),mrs)
               call mcstvt(hess(llhcc),xk(lxk),bijt(mk,mk),nbr,nbs)
c
 490        continue
            icc=icc+1
            llhcc=llhcc+mrs
            if(icc.eq.nbufcc.and.incorh.eq.0) then
c-2        write(iout,*)'  write hess_cc  hess ',nwcc,hess(lhcc)
               llhcc=lhcc
               call iosys('write real hess_cc on rwf without rewinding',
     $              nwcc,hess(llhcc),0,' ')
               call rzero(hess(llhcc),nwcc)
               icc=0
            endif

            lxk=lxk+nrs
 496     lxj=lxj+nnp
c
 500  continue
c
      if(iac.ne.0.and.incorh.eq.0) then
c-2       write(iout,*)' final write to hess_ac '
         nlwac=iac*mrs
         call iosys('write real hess_ac on rwf without rewinding',
     $        nlwac,hess(lhac),0,' ')
      endif
c
      if(icc.ne.0.and.incorh.eq.0) then
c-2       write(iout,*)' final write to hess_cc '
         nlwcc=icc*mrs
         call iosys('write real hess_cc on rwf without rewinding',
     $        nlwcc,hess(lhcc),0,' ')
      endif
c
      return
      end
