*deck %W%  %G%
      subroutine natodr(cr,c,nwwp,nfo,nco,nao,nvo,nob,nbf,
     $     nsym,nato,eigval,ifock,opnscf,iobcas)
c
c***begin prologue     natodr
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
c***end prologue       natodr
c
c
      implicit real*8(a-h,o-z)
c
      dimension cr(2),c(2)
c...
c...  real*8 cr(nwwp)
c...
c
      dimension nob(10),nco(10),nfo(10),nao(10),nvo(10),
     $     nbf(10),nato(10),lfab(10),iobcas(*)
c
      integer opnscf,wpadti
c
c
      ivt=0
      nbft=0
      ldent=0
      ixx=1
      maxnbf=0
      do 10 i=1,nsym
         lfab(i)=ixx
         nobi=nob(i)
         nbfi=nbf(i)
         nbf2=nbfi*nbfi
         ldent=ldent+(nao(i)*(nao(i)+1))/2
         ixx=ixx+nbf2
         nbft=nbft+nbfi
         maxnbf=max(nbf2,maxnbf)
         ivt=ivt+nbf2
 10   continue
c
c     allocate core
      ic=1
      if=ic+ivt
      ict=if+ivt+2
      iden=ict+ivt
      itemp=iden+ldent
      itv=itemp+maxnbf
      ievec=itv+maxnbf
      ieval=ievec+maxnbf
      ineed=wpadti(ieval+nbft)
c..bhl..unicos
      call getscm(ineed,cr,ngot,'natodr',0)
cc
      icc=1
      do 99 i=1,nsym
         ntotm=nob(i)*nbf(i)
         call vmove(cr(ict-1+icc),c(icc),ntotm)
         icc=icc+nbf(i)*nob(i)
 99   continue
cc
      call natral(lfab,cr(if),cr(iden),cr(itemp),cr(itv),
     $     cr(ict),c,cr(ievec),eigval,
     $     nsym,nfo,nco,nao,nvo,nbf,nato,ifock,opnscf,iobcas)
cc
      return
      end
