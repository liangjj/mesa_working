*deck @(#)trn2e.f	5.1  11/6/94
      subroutine trn2e(values,nnp,ntriang,c,nbf,nco,nocc,t1,t2,numij,
     $     valj,valk, nbufj,nbufk, val,t3,nao,lst16,kstrt)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      871116   (yymmdd)
c
c  16 november 1987    bhl at lanl
c   return inserted if nao = 0
c
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)trn2e.f	5.1   11/6/94
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
c
c***purpose: no-symmetry four index transformation of two-electron
c            integrals.
c
c paul saxe                 21 august 1984                  lanl
c
c bhl   array t1 and t3 are equivalenced in the call to this routine
c bhl   array val and values are equivalenced in a like manner
c
      implicit integer (a-z)
c
      real*8 values(nnp,ntriang),c(nbf,nocc),t1(nocc,nocc)
      real*8 t2(nocc,nbf),sqrt2,t3(nao,nao)
      real*8 val(numij,numij),valj(2),valk(2),cnk,cnl
      character*3 ians
c
      common /io/ inp,iout
c
c
      sqrt2=sqrt(2.0d+00)
      nbf22=nbf*nbf
c
c     ----- find where the first location in a word addressable file is -----
c
      call iosys('rewind "sorted ao integrals" on ints',0,0,0,' ')
c
c.io      call srew(itap91)
c.io      call getwad(itap91,st91)
c
c     ----- tape2 is positioned at the end of the label already -----
c
c      call getwad(itape2,st2)
c
c     ----- loop through kl triangles of integrals, transforming -----
c
cc
 51   format(//,'  first half transform ',/)
cc
      nbf22=nbf*nbf
      nnp=nbf*(nbf+1)/2
      nbf2=nnp
      nnb=(nocc*(nocc+1))/2
      lvalj=nnb*nnp
      lvalk=nnb*nbf22
c
      call iosys('does mc_j exist on mcscr',0,0,0,ians)
c
      if(ians.eq.'no') then
         mtotj=nnb*nbf2+2
         call iosys('create real mc_j on mcscr',mtotj,0,0,' ')
         mtotk=nnb*nbf22+2
         call iosys('create real mc_k on mcscr',mtotk,0,0,' ')
      endif
cc
c
      do 1 i=1,lvalj
         valj(i)=0.0d+00
 1    continue
      do 2 i=1,lvalk
         valk(i)=0.0d+00
 2    continue
cc
      maxkl=0
      kl=0
      klc=0
      do 4 k=1,nbf
         do 3 l=1,k
            kl=kl+1
            klc=klc+1
c
c     ----- check that this triangle of integrals is in core -----
c
            if (kl.gt.maxkl) then
               minkl=maxkl+1
               maxkl=min(nnp,maxkl+ntriang)
               lnread=(maxkl-minkl+1)*nnp
c               call sread(itap44,values,lnread)
               call iosys('read real "sorted ao integrals" from ints'//
     $              ' without rewinding',lnread,values,0,' ')
            end if
c
c     ----- perform first half-transformation -----
c
            call trtosq(t1,values(1,kl-minkl+1),nbf,nnp)
cc
cmp            write(iout,11119)
11119       format(/,' raw  ')
cmp         call prtint(t1,nbf,k,l)
cc
            call ebtc(t2,c,t1,nocc,nbf,nbf)
            call ebc (t1,t2,c,nocc,nbf,nocc)
c
c            write(iout,11118)
11118       format(/,' half   ')
c
c         call prtint(t1,nbf,nbf,l)
c
            ij=klc
            do 22 i=1,nocc
               do 23 j=1,i
                  valj(ij)=t1(i,j)
                  ij=ij+nnp
 23            continue
 22         continue
c-----------------------------------------------------------c
c
c     ---- exchange tranform  ----                          c
c
c       valk(m,n,j,l) = valk(m,n,j,l) + c(k,m)*t2(k,l,n,j)
c
c-----------------------------------------------------------c
            nx=0
            if(k.eq.l) go to 28
            nkk=(k-1)*nbf+1
            nll=(l-1)*nbf+1
            do 27 m=1,nocc
               do 26 n=1,m
                  cnl=c(l,n)
                  cnk=c(k,n)
                  jk=nkk
                  jl=nll
                  do 25 j=1,nbf
                     valk(jk)=valk(jk)+cnl*t2(m,j)
                     valk(jl)=valk(jl)+cnk*t2(m,j)
                     jl=jl+1
                     jk=jk+1
 25               continue
                  nkk=nkk+nbf22
                  nll=nll+nbf22
 26            continue
 27         continue
c
            go to 3
c
 28         continue
            nll=(l-1)*nbf+1
            do 31 m=1,nocc
               do 30 n=1,m
                  cnk=c(k,n)
                  jl=nll
                  do 29 j=1,nbf
                     valk(jl)=valk(jl)+cnk*t2(m,j)
                     jl=jl+1
 29               continue
                  nll=nll+nbf22
 30            continue
 31         continue
c
c
c
c
    3    continue
    4 continue
c
      lst16=0
      iblock=nbufj
c
         call iosys('write real mc_j on mcscr',
     $        lvalj,valj,lst16,' ')
c
      kstrt=1
      kst16=0
c
         call iosys('write real mc_k on mcscr',
     $        lvalk,valk,kst16,' ')
c
c
c     ----- loop through ij triangles of integrals, transforming -----
c
c     ----- only active-active block is transformed
c
c..bhl..closed shell scf
c
      if(nao.eq.0) return
c
      pt=0
      ij=0
      npt=1
      do 33 i=1,nao
         do 32 j=1,i
            ij=ij+1
c
            call trtosq(t3,valj(npt),nbf,nnp)
            npt=npt+nnp
c
            call ebtc(t2,c,t3,nao,nbf,nbf)
            call ebc (t3,t2,c,nao,nbf,nao)
cc
cc
            ijt=0
            do 40 it=1,nao
               do 41 jt=1,it
                  ijt=ijt+1
                  val(ijt,ij)=t3(it,jt)
 41            continue
 40         continue
cc
cc
 32      continue
 33   continue
c
c
      return
      end
