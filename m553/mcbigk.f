*deck @(#)mcbigk.f	5.1  11/6/94
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
c***source              @(#)mcbigk.f	5.1   11/6/94
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
      logical debug
c
      character*3 ians
c
      common /io/ inp,iout
      common / number / zero,pt5,one,two,four,eight
c
      parameter (debug=.false.)
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
      if(debug) then
         write(iout,*) 'before mclden'
      endif
      call mclden(ipqrs,mk,ml,nd,nstep,nij,mjk,
     $            ntot,ittype,nu,isqr,nblkd,mblkd,lblkd)
c
      if(ittype.ne.idtype) then
         write(iout,409)
 409     format('  mcdnk  idtype ne ittype .. stop')
         call lnkerr(' ')
      endif
c
c
      nleft=0
c
      ipr=isymc(mr)
      ips=isymc(ms)
      ixab=0
      if(noci.eq.0.and.nco(mr).ne.0.and.(nao(mr)+nvirt(mr)).ne.0)  then
         ixab=1
         fctiv1=-two*deg(mr)
         fctiv2=-one*deg(mr)
         n34=(nao(mk)*(nao(mk)+1))/2
         n12=(nco(mr)*(nco(mr)+1))/2
         ibnr=nco(1)*nob(1)
         ibnnk=(nao(1)+1)*nao(1)/2
         lpass=lbufso/ibnr
         npass=(ibnnk-1)/lpass+1
c
         call iosys('does abix exist on rwf',0,0,0,ians)
         if (ians.eq.'no') then
            call iosys('create real abix on rwf',npass*lbufso,0,0,' ')
         end if
         call iosys('rewind abix on mcscr',0,0,0,' ')
         call iosys('rewind abix on rwf',0,0,0,' ')
         call iosys('read real abix from mcscr without rewinding',
     $              lbufso,bufab,0,' ')
         kpass=1
c
      endif
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
            if(iread.gt.iijt) then
               ileft=min(iijt,noc2-ird)
               nbufkk=ileft*nrs
               ird=ird+ileft
               lxk=1
               iread=1
               kstsav=kst
               call mcrdk(xk,nbufkk,ndab,kst,jkcore)
            endif
c
c
            if(nstep.eq.1) then
c
c------------------------------------------------------
c   multiply the exchange operator times the appropiate
c   density matrix elements
c------------------------------------------------------
c
               id =((l-1)*naol+k-1)*nij
               if(debug) then
                  write(iout,*) 'before mcdxaa'
               endif
               call mcdxaa(den(id+nd),xk(lxk),nij,mrs,nstep,itran,thrsh,
     $                     hess(lhaa),nstart,ntot)
c
               if(k.ne.l) then
c
c---------------------------------------------------------
c   the transpose of this exchange operator is also needed
c---------------------------------------------------------
c
                  idt=((k-1)*naol+l-1)*nij
                  if(debug) then
                     write(iout,*) 'before mctrsp'
                  endif
                  call mctrsp(xk(lxk),xkt,nbr,nbs)
                  if(debug) then
                     write(iout,*) 'before mcdxaa'
                  endif
                  call mcdxaa(den(idt+nd),xkt,nij,mrs,nstep,itran,thrsh,
     $                        hess(lhaa),nstrtt,ntott)
cc
               endif
            else
               write(iout,416)
 416           format(' problem mcdnk transpose d.m. ? ')
               call lnkerr(' ')
            endif
c
c-------------------------------------------
c   xk contributes to the fab  fock-operator
c-------------------------------------------
c
            ldab1=ldab(mk)+(k-1)*naol+l-1
            lfab1=lfab(mr)
            if(k.ne.l) then
c
               factk=-deg(mr)*dab(ldab1)
c
               if(debug) then
                  write(iout,*) 'before mcstv'
               endif
               call mcstv(fab(lfab1),xk(lxk),factk,mrs)
            else
               factk=-pt5*deg(mr)*dab(ldab1)
               if(debug) then
                  write(iout,*) 'before mcstv'
               endif
               call mcstv(fab(lfab1), xk(lxk),factk,mrs)
            endif
c
c--------------------------------------
c   core contribution to abii  i-vector
c--------------------------------------
c
            if(ixab.ne.0) then
               if(k.ne.l) then
                  ikt=1
                  if(debug) then
                     write(iout,*) 'before mcaibx'
                  endif
                  call mcaibx(c(ipr),xk(lxk),nbr,nco(mr),nob(mr),temp,
     $                       tv,nfiv,fctiv2,ikt,'aibx',bufix(ibix+1),
     $                       bufab(ibix+1))
               else
                  ikt=0
                  if(debug) then
                     write(iout,*) 'before mcaibx'
                  endif
                  call mcaibx(c(ipr),xk(lxk),nbr,nco(mr),nob(mr),
     $                        temp,tv,nfiv,fctiv1,ikt,'aibx',
     $                        bufix(ibix+1),bufab(ibix+1))
               endif
               ibix=ibix+lenix
               if (ibix+lenix.gt.lbufso) then
                  call iosys('write real abix to rwf without rewinding',
     $                       lbufso,bufab,0,' ')
                  if(kpass.lt.npass) then
                     call iosys('read real abix from mcscr '//
     $                          'without rewinding',lbufso,bufab,0,' ')
                     kpass=kpass+1
                  end if
                  ibix=0
               end if
c
c
            endif
            lxk=lxk+nrs
  495    continue
 500  continue
c
c
      if (ibix.ne.0) then
         call iosys('write real abix to rwf without rewinding',
     $              lbufso,bufab,0,' ')
      end if
c
c
      return
      end
