*deck @(#)xxint.f	5.1  11/6/94
      subroutine xxint(int,h,orbsym,ijgrp,ijadd,kadd,ladd,ijww,klww,
     #                 ijxx,klxx,val,lab,bin,offset,minsym,maxsym,
     #                 lenbin,asort,prnt)
c
      implicit real*8 (a-h,o-z)
c
      integer xor
      logical prnt
      integer symorb
      integer maxb,lvfrm1,pt
      real*8 int(nmax),h(nijvir),val(lenbin),asort(*)
      integer kadd(symorb),ladd(symorb),ijadd(numij),orbsym(norbs)
      integer ijgrp(numij),ijww(numij),klww(nijvir),ijxx(numij)
      integer klxx(nijvir),lab(lenbin),bin(lenbin),offset(norbs)
      integer minsym(0:nsym-1),maxsym(0:nsym-1)
c
      common /dimn/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,nrow4o,nwksoc,nlevoc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,maxb,nroots,lvfrm1,nrefs
      common /intm/ nmax,ngroup,nblkoc,numij,symorb,intsrt
      common /io/     inp,iout
      common /sort/   iblock,lnbuf,maxsrt
      common /x4x821/   nijvir
c
      sqrt2=sqrt(2.0d+00)
c
      do 1 i=1,norbs
         orbsym(i)=orbsym(i)-1
         offset(i)=i*(i-1)/2
    1 continue
      do 131 i=0,nsym-1
         minsym(i)=0
         maxsym(i)=-1
  131 continue
      is=orbsym(1)
      minsym(is)=1
      do 132 i=2,lvfrm1
         if (orbsym(i).eq.is) go to 132
            maxsym(is)=i-1
            is=orbsym(i)
            minsym(is)=i
  132 continue
      maxsym(is)=lvfrm1
c
c     ----- get the one-electron mo integrals -----
c
      call iosys('rewind "guga integrals" on gints',0,0,0,' ')
      call iosys('read real "guga integrals" from gints',nmax,int,0,' ')

      iblock=1
c
      do 6 i=lvfrm1,2,-1
         ia=i*(i-1)/2
         ii=ia+i
    2    if (ijgrp(ii).eq.iblock) go to 3
            iblock=iblock+1
            if (iblock.gt.ngroup) stop
            call iosys('read real "guga integrals" from gints '//
     #                 'without rewinding',nmax,int,0,' ')
         go to 2
c
    3    continue
         iii=ijadd(ia+i)+kadd(i)
         iiis=orbsym(i)
         iiin=iiis*norbs
         iiij=iii+ladd(iiin+minsym(iiis))+3
         do 5 j=minsym(iiis),min(i-1,maxsym(iiis))
            h(ia+j)=int(iiij)
            iiij=iiij+3
    5    continue
    6 continue
c
c
c
      call iosys('rewind "guga integrals" on gints',0,0,0,' ')
      call iosys('read real "guga integrals" from gints',nmax,int,0,' ')
      iblock=1
c
c     ----- form the 3-external supermatrices ------
c
      pt=0
      do 110 i=norbs,levfrm,-1
         ia=offset(i)
         is=orbsym(i)
         do 109 j=lvfrm1,2,-1
            ij=ia+j
  101       continue
               if (ijgrp(ij).eq.iblock) go to 102
               iblock=iblock+1
               if (iblock.gt.ngroup) stop 1
               call iosys('read real "guga integrals" from gints '//
     #                    'without rewinding',nmax,int,0,' ')
               go to 101
  102       continue
            ja=offset(j)
            ijs=xor(is,orbsym(j))
            ijn=ijs*norbs
            ijp=ijadd(ij)
c
c     ----- [ij;jl] and [il;jj] -----
c
            jj=ja+j
            ijj=ijp+kadd(ijn+j)
            ijjn=is*norbs
            minl=minsym(is)
            maxl=min(j-1,maxsym(is))
            ijjl=ijj+ladd(ijjn+minl)
            il=ia+minl
            jl=ja+minl
c
            if (pt+3*(maxl-minl+1).gt.lenbin) then
c              call sorter('with bin',asort,asort,0,0,0,0,
c    #                      0,0,0,0,val,lab,bin,pt,prnt)
            call sorter('with bin',asort,asort,0,pt,lab,bin,val,0,0,0,
     $                   prnt)
               pt=0
            end if
c
cdir$ ivdep
            do 104 l=minl,maxl
               t1=int(ijjl+1)
               t2=int(ijjl+2)
               val(pt+1)=t1*sqrt2
               lab(pt+1)=ijww(il)+klww(jj)
               val(pt+2)=t1+t2
               lab(pt+2)=ijww(ij)+klww(jl)
               val(pt+3)=t1-t2
               lab(pt+3)=ijxx(ij)+klxx(jl)
               pt=pt+3
               ijjl=ijjl+3
               il=il+1
               jl=jl+1
  104       continue
c
            do 108 k=1,j-1
               ijk=ijp+kadd(ijn+k)
               ik=ia+k
               jk=ja+k
               ka=offset(k)
               ijks=xor(ijs,orbsym(k))
               ijkn=ijks*norbs
c
c     ----- [ij;kk] and [ik;jk] -----
c
               if (orbsym(k).eq.ijks) then
                  ijkk=ijk+ladd(ijkn+k)
                  kk=ka+k
c
                  if (pt+3.gt.lenbin) then
c                    call sorter('with bin',asort,asort,0,0,0,0,
c    #                            0,0,0,0,val,lab,bin,pt,prnt)
            call sorter('with bin',asort,asort,0,pt,lab,bin,val,0,0,0,
     $                   prnt)
                     pt=0
                  end if
c
                  t1=int(ijkk+1)
                  t2=int(ijkk+2)
                  val(pt+1)=t1*sqrt2
                  lab(pt+1)=ijww(ij)+klww(kk)
                  val(pt+2)=t2+t1
                  lab(pt+2)=ijww(ik)+klww(jk)
                  val(pt+3)=t2-t1
                  lab(pt+3)=ijxx(ik)+klxx(jk)
                  pt=pt+3
               end if
c
c     ----- [ij;kl] and [ik;jl] and [il;jk] -----
c
               minl=minsym(ijks)
               maxl=min(k-1,maxsym(ijks))
               ijkl=ijk+ladd(ijkn+minl)
               il=ia+minl
               jl=ja+minl
               kl=ka+minl
c
               if (pt+6*(maxl-minl+1).gt.lenbin) then
c                 call sorter('with bin',asort,asort,0,0,0,0,
c    #                         0,0,0,0,val,lab,bin,pt,prnt)
            call sorter('with bin',asort,asort,0,pt,lab,bin,val,0,0,0,
     $                   prnt)
                  pt=0
               end if
c
cdir$ ivdep
               do 107 l=minl,maxl
                  t1=int(ijkl+1)
                  t2=int(ijkl+2)
                  t3=int(ijkl+3)
                  val(pt+1)=t1+t3
                  val(pt+2)=t1-t3
                  val(pt+3)=t2+t3
                  val(pt+4)=t2-t3
                  val(pt+5)=t2+t1
                  val(pt+6)=t2-t1
                  lab(pt+1)=ijww(ij)+klww(kl)
                  lab(pt+2)=ijxx(ij)+klxx(kl)
                  lab(pt+3)=ijww(ik)+klww(jl)
                  lab(pt+4)=ijxx(ik)+klxx(jl)
                  lab(pt+5)=ijww(il)+klww(jk)
                  lab(pt+6)=ijxx(il)+klxx(jk)
                  pt=pt+6
                  ijkl=ijkl+3
                  il=il+1
                  jl=jl+1
                  kl=kl+1
  107          continue
  108       continue
  109    continue
  110 continue
c
c     ----- four-external integrals -----
c
      do 23 i=lvfrm1,2,-1
         ia=offset(i)
         is=orbsym(i)
         ii=ia+i
c
    7    continue
            if (ijgrp(ii).eq.iblock) go to 8
            iblock=iblock+1
            if (iblock.gt.ngroup) stop 1
            call iosys('read real "guga integrals" from gints '//
     #                 'without rewinding',nmax,int,0,' ')
            go to 7
c
    8    continue
         iip=ijadd(ii)
         iii=iip+kadd(i)
         iiin=is*norbs
c
c     ----- [ii;il] and [il;ll] -----
c
         minl=minsym(is)
         maxl=min(i-1,maxsym(is))
         il=ia+minl
         iiil=iii+ladd(iiin+minl)
c
         if (pt+4*(maxl-minl+1).gt.lenbin) then
c           call sorter('with bin',asort,asort,0,0,0,0,
c    #                   0,0,0,0,val,lab,bin,pt,prnt)
            call sorter('with bin',asort,asort,0,pt,lab,bin,val,0,0,0,
     $                   prnt)
            pt=0
         end if
c
cdir$ ivdep
         do 10 l=minl,maxl
            ll=offset(l)+l
            val(pt+1)=sqrt2*(int(iiil+1)+int(iiil+3))
            val(pt+2)=val(pt+1)
            lab(pt+1)=ijww(ii)+klww(il)
            lab(pt+2)=ijww(il)+klww(ii)
            val(pt+3)=sqrt2*(int(iiil+2)+int(iiil+3))
            val(pt+4)=val(pt+3)
            lab(pt+3)=ijww(il)+klww(ll)
            lab(pt+4)=ijww(ll)+klww(il)
            pt=pt+4
            il=il+1
            iiil=iiil+3
   10    continue
c
c     ----- [ii;kk] -----
c
         do 13 k=1,i-1
            iik=iip+kadd(k)
            iiks=orbsym(k)
            iikn=iiks*norbs
            ka=offset(k)
c
            ik=ia+k
            kk=ka+k
            iikk=iik+ladd(iikn+k)
c
            if (pt+2.gt.lenbin) then
c              call sorter('with bin',asort,asort,0,0,0,0,
c    #                      0,0,0,0,val,lab,bin,pt,prnt)
            call sorter('with bin',asort,asort,0,pt,lab,bin,val,0,0,0,
     $                   prnt)
               pt=0
            end if
c
            val(pt+1)=int(iikk+1)
            val(pt+2)=int(iikk+1)
            lab(pt+1)=ijww(ii)+klww(kk)
            lab(pt+2)=ijww(kk)+klww(ii)
            pt=pt+2
c
c     ----- [ii;kl] and [ik;il] -----
c
            minl=minsym(iiks)
            maxl=min(k-1,maxsym(iiks))
            il=ia+minl
            kl=ka+minl
            iikl=iik+ladd(iikn+minl)
c
            if (pt+6*(maxl-minl+1).gt.lenbin) then
c              call sorter('with bin',asort,asort,0,0,0,0,
c    #                      0,0,0,0,val,lab,bin,pt,prnt)
            call sorter('with bin',asort,asort,0,pt,lab,bin,val,0,0,0,
     $                   prnt)
               pt=0
            end if
c
cdir$ ivdep
            do 12 l=minl,maxl
               t1=int(iikl+1)
               t2=int(iikl+2)
               hkl=h(kl)
               val(pt+1)=t2-t1+hkl
               val(pt+2)=val(pt+1)
               lab(pt+1)=ijxx(ik)+klxx(il)
               lab(pt+2)=ijxx(il)+klxx(ik)
               val(pt+3)=sqrt2*t1
               val(pt+4)=val(pt+3)
               lab(pt+3)=ijww(ii)+klww(kl)
               lab(pt+4)=ijww(kl)+klww(ii)
               val(pt+5)=t1+t2+hkl
               val(pt+6)=val(pt+5)
               lab(pt+5)=ijww(ik)+klww(il)
               lab(pt+6)=ijww(il)+klww(ik)
               pt=pt+6
               il=il+1
               kl=kl+1
               iikl=iikl+3
   12       continue
   13    continue
c
         do 22 j=i-1,2,-1
            ij=ia+j
c
   14       continue
               if (ijgrp(ij).eq.iblock) go to 15
               iblock=iblock+1
               if (iblock.gt.ngroup) stop 1
      call iosys('read real "guga integrals" from gints '
     $                    //'without rewinding',nmax,int,0,' ')
               go to 14
c
   15       continue
            hij=h(ij)
            ijp=ijadd(ij)
            ijs=xor(is,orbsym(j))
            ijn=ijs*norbs
            ja=offset(j)
            jj=ja+j
c
c     ----- [ij;jl] and [il;jj] -----
c
            ijj=ijp+kadd(ijn+j)
            ijjs=is
            ijjn=ijjs*norbs
            minl=minsym(ijjs)
            maxl=min(j-1,maxsym(ijjs))
            il=ia+minl
            jl=ja+minl
            ijjl=ijj+ladd(ijjn+minl)
c
            if (pt+6*(maxl-minl+1).gt.lenbin) then
c              call sorter('with bin',asort,asort,0,0,0,0,
c    #                      0,0,0,0,val,lab,bin,pt,prnt)
            call sorter('with bin',asort,asort,0,pt,lab,bin,val,0,0,0,
     $                   prnt)
               pt=0
            end if
c
cdir$ ivdep
            do 17 l=minl,maxl
               t1=int(ijjl+1)
               t2=int(ijjl+2)
               hil=h(il)
               val(pt+1)=t1-t2-hil
               val(pt+2)=val(pt+1)
               lab(pt+1)=ijxx(ij)+klxx(jl)
               lab(pt+2)=ijxx(jl)+klxx(ij)
               val(pt+3)=sqrt2*t1
               val(pt+4)=val(pt+3)
               lab(pt+3)=ijww(il)+klww(jj)
               lab(pt+4)=ijww(jj)+klww(il)
               val(pt+5)=t1+t2+hil
               val(pt+6)=val(pt+5)
               lab(pt+5)=ijww(ij)+klww(jl)
               lab(pt+6)=ijww(jl)+klww(ij)
               pt=pt+6
               il=il+1
               jl=jl+1
               ijjl=ijjl+3
   17       continue
c
c     ----- [ij;kk] and [ik;jk] -----
c
            do 21 k=1,j-1
               ijk=ijp+kadd(ijn+k)
               ijks=xor(ijs,orbsym(k))
               ijkn=ijks*norbs
               ik=ia+k
               jk=ja+k
               ka=offset(k)
c
               if (orbsym(k).eq.ijks) then
                  kk=ka+k
                  ijkk=ijk+ladd(ijkn+k)
                  if (pt+6.gt.lenbin) then
c                    call sorter('with bin',asort,asort,0,0,0,0,
c    #                            0,0,0,0,val,lab,bin,pt,prnt)
            call sorter('with bin',asort,asort,0,pt,lab,bin,val,0,0,0,
     $                   prnt)
                     pt=0
                  end if
                  t1=int(ijkk+1)
                  t2=int(ijkk+2)
                  val(pt+1)=t2-t1+hij
                  val(pt+2)=val(pt+1)
                  lab(pt+1)=ijxx(ik)+klxx(jk)
                  lab(pt+2)=ijxx(jk)+klxx(ik)
                  val(pt+3)=sqrt2*t1
                  val(pt+4)=val(pt+3)
                  lab(pt+3)=ijww(ij)+klww(kk)
                  lab(pt+4)=ijww(kk)+klww(ij)
                  val(pt+5)=t1+t2+hij
                  val(pt+6)=val(pt+5)
                  lab(pt+5)=ijww(ik)+klww(jk)
                  lab(pt+6)=ijww(jk)+klww(ik)
                  pt=pt+6
               end if
c
c     ----- [ij;kl] and [ik;jl] and [il;jk] -----
c
               minl=minsym(ijks)
               maxl=min(k-1,maxsym(ijks))
               il=ia+minl
               jl=ja+minl
               kl=ka+minl
               ijkl=ijk+ladd(ijkn+minl)
c
               if (pt+12*(maxl-minl+1).gt.lenbin) then
c                 call sorter('with bin',asort,asort,0,0,0,0,
c    #                         0,0,0,0,val,lab,bin,pt,prnt)
            call sorter('with bin',asort,asort,0,pt,lab,bin,val,0,0,0,
     $                   prnt)
                  pt=0
               end if
c
cdir$ ivdep
               do 20 l=minl,maxl
                  t1=int(ijkl+1)
                  t2=int(ijkl+2)
                  t3=int(ijkl+3)
                  val(pt+1)=t1-t3
                  val(pt+2)=val(pt+1)
                  lab(pt+1)=ijxx(ij)+klxx(kl)
                  lab(pt+2)=ijxx(kl)+klxx(ij)
                  val(pt+3)=t1+t3
                  val(pt+4)=val(pt+3)
                  lab(pt+3)=ijww(ij)+klww(kl)
                  lab(pt+4)=ijww(kl)+klww(ij)
                  val(pt+5)=t2-t3
                  val(pt+6)=val(pt+5)
                  lab(pt+5)=ijxx(ik)+klxx(jl)
                  lab(pt+6)=ijxx(jl)+klxx(ik)
                  val(pt+7)=t2+t3
                  val(pt+8)=val(pt+7)
                  lab(pt+7)=ijww(ik)+klww(jl)
                  lab(pt+8)=ijww(jl)+klww(ik)
                  val(pt+9)=t2-t1
                  val(pt+10)=val(pt+9)
                  lab(pt+9)=ijxx(il)+klxx(jk)
                  lab(pt+10)=ijxx(jk)+klxx(il)
                  val(pt+11)=t2+t1
                  val(pt+12)=val(pt+11)
                  lab(pt+11)=ijww(il)+klww(jk)
                  lab(pt+12)=ijww(jk)+klww(il)
                  pt=pt+12
                  il=il+1
                  jl=jl+1
                  kl=kl+1
                  ijkl=ijkl+3
   20          continue
   21       continue
   22    continue
   23 continue
c
c     ----- flush the sort bins -----
c
c     call sorter('with bin',asort,asort,0,0,0,0,
c    #             0,0,0,0,val,lab,bin,pt,prnt)
            call sorter('with bin',asort,asort,0,pt,lab,bin,val,0,0,0,
     $                   prnt)
c
c
      return
      end
