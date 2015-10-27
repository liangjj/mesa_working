*deck @(#)jtran.f	5.1  11/6/94
      subroutine jtran(values,nnp,ntriang,c,nbf,nco,nao,nocc,
     $     t1,t2,t3,
     $     asort,ncors, valt,xval, lsortj,
     $     xjbuf,jlab,jjlab,labj,nbufj,
     $     xkbuf,klab,kklab,labk,nbufk)

c
c***begin prologue     jtran
c***date written       871120   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (lanl)
c***source             @(#)jtran.f	5.1   11/6/94
c
c***purpose            no-symmetry two index transformation of
c                      two-electron integrals to form coulomb
c                      and exchange intregral block. the exchange
c                      transformation is completed in routine ktran
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       jtran
c
      implicit integer (a-z)
c
      real*8 values(nnp,ntriang),c(nbf,nocc),t1(2),valt(2),xval(2)
      real*8 t2(2),sqrt2,t3(2),asort(*)
      real*8 xjbuf(nbufj),xkbuf(nbufk)
      dimension jlab(nbufj),jjlab(nbufj),labj(*),klab(nbufk)
      integer kklab(nbufk),labk(*)

c
      common /io/ inp,iout
c
      sqrt2=sqrt(2.0d+00)
c
c     ----- find where the first location in a word addressable
c                                                         file is -----
c
      call iosys('rewind "sorted ao integrals" on ints',0,0,0,' ')
c
c     ----- loop through kl triangles of integrals, transforming -----
c
      nbf22=nbf*nbf
      nbf2=nbf*(nbf+1)/2
      nocnbf=nocc*nbf
      noc2=(nocc*(nocc+1))/2
      nna=(nao*(nao+1))/2
      lvalt=nna*nao*nbf
      call rzero(valt,lvalt)
c
      lsortj=noc2*nbf2
      lsortk=nocnbf*nbf2
      call sorter('start',asort,asort,wpadti(ncors),lsortk+lsortj+1,
     $     0,0,0,0,'mc_j','mcscr',.false.)
c
      no=0
      do 1 i=1,noc2
      labj(i)=no
      no=no+nbf2
 1    continue
c
      no=lsortj
      do 2 i=1,nocnbf
      labk(i)=no
      no=no+nbf2
 2    continue
c
c     loop over passes
c
      maxkl=0
      kl=0
      do 1000 k=1,nbf
         do 999  l=1,k
c
            kl=kl+1
c
c     ----- check that this triangle of integrals is in core -----
c
            if (kl.gt.maxkl) then
               minkl=maxkl+1
               maxkl=min(nnp,maxkl+ntriang)
               lnread=(maxkl-minkl+1)*nnp
               call iosys('read real "sorted ao integrals" from ints'//
     $              ' without rewinding',lnread,values,0,' ')
            end if
c
c     ----- perform first half-transformation -----
c
            call trtosq(t1,values(1,kl-minkl+1),nbf,nnp)
            call ebtc(xkbuf,c,t1,nocc,nbf,nbf)
            call ebc (t1,xkbuf,c,nocc,nbf,nocc)
            call trabcx(c,valt,t1,nbf,nco,nao,nocc,k,l)
            call sqtotr(xjbuf,t1,nocc,noc2)
cccccccccc
c
c     ---- labels for coulomb sort ----
c
cccccccccc
            do 50 n=1,noc2
               jlab(n)=labj(n)+kl
 50         continue
cccccccccc
c
c     ---- labels for exchange sort ----
c
cccccccccc
            do 60 n=1,nocnbf
               klab(n)=labk(n)+kl
 60         continue
c
c     ---- sort the ints  ----
c
            call sorter('with bin',asort,asort,0,noc2,jlab,jjlab,
     $           xjbuf,0,0,0,.false.)
c
            call sorter('with bin',asort,asort,0,nocnbf,klab,kklab,
     $           xkbuf,0,0,0,.false.)
c
 999     continue
 1000 continue
c
c ---- finish the  sort
c
          call sorter('end',asort,asort,0,0,0,0,0,0,0,0,.false.)
c
c ---- finish the active-active-active-active transform
c
      nao3=nao*nna
      call ebtc(xval,c,valt,nao,nbf,nao3)
      call fixval(values,xval,nao)
c
      return
      end
