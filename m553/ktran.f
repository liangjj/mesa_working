*deck @(#)ktran.f	5.1  11/6/94
      subroutine ktran(values,nnp,ntriang, c, nbf,nco,nao,nocc, t1,
     $     asort,ncors,  lsortj,
     $     xkbuf,klab,kklab,labij,nbufk)

c
c***begin prologue     ktran
c***date written       871120   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (lanl)
c***source             @(#)ktran.f	5.1   11/6/94
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
c***end prologue       ktran
c
      implicit integer (a-z)
c
      real*8 values(nnp,ntriang),c(nbf,nocc),t1(*)
      real*8 xkbuf(nbufk),asort(*)
      dimension labij(nocc,nocc),klab(nbufk),kklab(nbufk)

c
      common /io/ inp,iout
c
c     ----- loop through kl triangles of integrals, transforming -----
c
      nbf22=nbf*nbf
      nbf2=nbf*(nbf+1)/2
      nocnbf=nocc*nbf
      noc2=(nocc*(nocc+1))/2
      nna=(nao*(nao+1))/2
c
      lsortk=noc2*nbf22
      call sorter('start',asort,asort,wpadti(ncors),lsortk+1,0,0,0,0,
     $ 'mc_k','mcscr',.false.)
c
      ix=0
      do 1 i=1,nocc
      do 2 j=1,i
      labij(i,j)=((ix+j)-1)*nbf22
 2    continue
      ix=ix+i
 1    continue
c
c     loop over passes
c
      maxkl=0
      kl=0
c
      kword=lsortj
c
      do 1000 k=1,nbf
         do 999  l=1,nocc
c
            kl=kl+1
c
c     ----- check that this triangle of integrals is in core -----
c
            if (kl.gt.maxkl) then
               minkl=maxkl+1
               maxkl=min(nocnbf,maxkl+ntriang)
               lnread=(maxkl-minkl+1)*nnp
               call iosys('read real mc_j from mcscr',
     $              lnread,values,kword,' ')
               kword=kword+lnread
            end if
c
c     ----- perform first half-transformation -----
c
            call trtosq(t1,values(1,kl-minkl+1),nbf,nnp)
            call ebc(xkbuf,t1,c,nbf,nbf,l)
cccccccccc
c
c     ---- labels for final exchange  sort ----
c
cccccccccc
            ix=0
            do 50 n=1,l
                labn=labij(l,n)+k
             do 60 m=1,nbf
               klab(ix+m)=labn
               labn=labn+nbf
 60          continue
             ix=ix+nbf
 50         continue
c
c     ---- sort the ints  ----
c
            call sorter('with bin',asort,asort,0,l*nbf,klab,kklab,
     $           xkbuf,0,0,0,.false.)
c
 999     continue
 1000 continue
c
c ---- finish the  sort
c
          call sorter('end',asort,asort,0,0,0,0,0,0,0,0,.false.)
c
      return
      end
