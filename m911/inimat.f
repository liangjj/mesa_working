*deck @(#)inimat.f	5.1  11/6/94
      subroutine inimat(int,ijadd,kadd,ladd)
c..rlm  subroutine inimat(int,ijadd,kadd,ladd,hiim,hiism,hijm,him,hsm)
c
c
c
      implicit real*8 (a-h,o-z)
c
      integer symorb,arr,orbfrm,asm,os,bmax
c
      real*8 int(nmax),ijadd(numij),kadd(symorb),ladd(symorb)
      real*8 hiim(1),hiism(1),hijm(1),him(1),hsm(1)
      common /coefs/  a,b,intal,intau,intb,intad,intbd
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
      common /dims/   nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     #,               nrowoc,nrow4o,nwksoc,nlevoc,norboc,levfrm
     #,               nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /all/    val1,val2,val3,arr,itr1,itr2,ia,ja,itype,isegt
     #,               lvfrm1,nlwki,nlwkj,imax,imin
      common /symm/   asm,js,is,mx(8),mn(8),os(8),numsym(8)
c
c
      data acrcy /1.0d-09/
      save acrcy
c
      return
c
c*********************************** hii *******************************
c
      entry hii(hiim,isym,numi,int,kadd)
c
      lad=arr+kadd(os(asm+1)+mn(isym+1))
c
c$    if (abs(a).lt.acrcy) go to 110
c$    if (abs(b).lt.acrcy) go to 120
c
      jisv=1
      do 102 i=1,numi
         ji=jisv
         ij=i
         do 101 j=1,i-1
            int(lad+intb )=int(lad+intb )+b*(hiim(ij)+hiim(ji))
            int(lad+intal)=int(lad+intal)+a*hiim(ij)
            int(lad+intau)=int(lad+intau)+a*hiim(ji)
            lad=lad+3
            ij=ij+numi
            ji=ji+1
  101    continue
         int(lad+intad)=int(lad+intad)+a*hiim(ij)
         int(lad+intbd)=int(lad+intbd)+b*hiim(ij)
         lad=lad+3
         jisv=jisv+numi
  102 continue
      return
cbl  dead code
c 110 continue
c     jisv=1
c     do 112 i=1,numi
c     ji=jisv
c     ij=i
c     do 111 j=1,i-1
c     t=b*int(lad+intb)
c     hiim(ij)=t
c     hiim(ji)=t
c     ji=ji+1
c     ij=ij+numi
c     lad=lad+3
c 111 continue
c     hiim(ij)=b*int(lad+intbd)
c     jisv=jisv+numi
c     lad=lad+3
c 112 continue
c     return
cblc dead code
c 120 continue
c     jisv=1
c     do 122 i=1,numi
c     ji=jisv
c     ij=i
c     do 121 j=1,i-1
c     hiim(ij)=a*int(lad+intal)
c     hiim(ji)=a*int(lad+intau)
c     ji=ji+1
c     ij=ij+numi
c     lad=lad+3
c 121 continue
c     hiim(ij)=a*int(lad+intad)
c     jisv=jisv+numi
c     lad=lad+3
c 122 continue
cc
c     return
c
c************************************ hiis *****************************
c
      entry hiis(hiism,isym,numi,int,kadd)
c
      lad=arr+kadd(os(asm+1)+mn(isym+1))
c
c$    if (abs(a).lt.acrcy) go to 210
c$    if (abs(b).lt.acrcy) go to 220
c
      jisv=1
      do 202 i=1,numi
         ji=jisv
         ij=i
         do 201 j=1,i-1
            t=hiism(ij)+hiism(ji)
            int(lad+intb )=int(lad+intb )+b*t
            int(lad+intal)=int(lad+intal)+a*t
            ji=ji+1
            ij=ij+numi
            lad=lad+3
  201    continue
         jisv=jisv+numi
         lad=lad+3
  202 continue
c
      return
cdeadcode
c 210 continue
c     jisv=1
c     do 212 i=1,numi
c     ji=jisv
c     ij=i
c     do 211 j=1,i-1
c     t=b*int(lad+intb)
c     hiism(ij)=t
c     hiism(ji)=t
c     ji=ji+1
c     ij=ij+numi
c     lad=lad+3
c 211 continue
c     hiism(ij)=0.0
c     jisv=jisv+numi
c     lad=lad+3
c 212 continue
c
c     return
c
c 220 continue
c     jisv=1
c     do 222 i=1,numi
c     ji=jisv
c     ij=i
c     do 221 j=1,i-1
c     t=a*int(lad+intal)
c     hiism(ij)=t
c     hiism(ji)=t
c     ji=ji+1
c     ij=ij+numi
c     lad=lad+3
c 221 continue
c     hiism(ij)=0.0
c     jisv=jisv+numi
c     lad=lad+3
c 222 continue
c
c     return
c
c********************************** hij ********************************
c
      entry hij(hijm,isym,jsym,numi,numj,int,kadd)
c
      lad=arr+kadd(os(asm+1)+mn(isym+1))
c
c$    if (abs(a).lt.acrcy) go to 310
c$    if (abs(b).lt.acrcy) go to 320
c
      do 302 i=1,numi
         ij=i
         do 301 j=1,numj
c           hijm(ij)=a*int(lad+intal)+b*int(lad+intb)
            int(lad+intal)=int(lad+intal)+a*hijm(ij)
            int(lad+intb )=int(lad+intb )+b*hijm(ij)
c
            ij=ij+numi
            lad=lad+3
  301    continue
  302 continue
c
      return
cdeadcode
c 310 continue
c     lad=lad+intb
c     do 312 i=1,numi
c     ij=i
c      do 311 j=1,numj
c     hijm(ij)=b*int(lad)
c     ij=ij+numi
c     lad=lad+3
c 311 continue
c 312 continue
c
c     return
c
c 320 continue
c     lad=lad+intal
c     do 322 i=1,numi
c     ij=i
c     do 321 j=1,numj
c     hijm(ij)=a*int(lad)
c     ij=ij+numi
c     lad=lad+3
c 321 continue
c 322 continue
c     return
c
c******************************* hi ************************************
c
      entry hi(him,isym,numi,int,ladd)
c
      lad=arr+ladd(os(asm+1)+mn(isym+1))
      do 400 i=1,numi
c        him(i)=a*int(lad+intal)+b*int(lad+intb)
         int(lad+intal)=int(lad+intal)+a*him(i)
         int(lad+intb )=int(lad+intb )+b*him(i)
         lad=lad+3
  400 continue
      return
c
c******************************* hs ************************************
c
      entry hs(hsm,isym,numi,int,ladd)
c
      a2=a*sqrt(2.0d+00)
      lad=arr+ladd(os(asm+1)+mn(isym+1))+2
      do 410 i=1,numi
c        hsm(i)=a2*int(lad)
         int(lad)=int(lad)+a2*hsm(i)
c
         lad=lad+3
  410 continue
      return
c
c
      end
