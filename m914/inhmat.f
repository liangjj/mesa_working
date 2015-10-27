*deck @(#)inhmat.f	1.1  11/30/90
      subroutine inhmat(int,ijadd,kadd,ladd,hiim,hiism,hijm,him,hsm)
csel  subroutine inimat(int,ijadd,kadd,ladd)
c
c
c
      implicit real*8 (a-h,o-z)
c
cvax  extended dummy int,ijadd,kadd,ladd,hiim,hiism,hijm,him,hsm
c
      integer symorb,arr,orbfrm,asm,os,bmax
      real*8 int
c
      common /coefs/  a,b,intal,intau,intb,intad,intbd
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
      common /dims/   nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     #,               nrowoc,levfrm
     #,               nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /all/    val1,val2,val3,arr,itr1,itr2,ia,ja,itype,isegt
     #,               lvfrm1,nlwki,nlwkj,imax,imin
      common /symq/   asm,js,is,mx(8),mn(8),os(8),numsym(8)
c
      dimension int(nmax),ijadd(numij),kadd(symorb),ladd(symorb)
      dimension hiim (1),hiism(1),hijm (1),him(1),hsm(1)
c
      data acrcy /1.0d-09/
c
      return
c
c************************************ inimat ***************************
c
      entry inimat(int,ijadd,kadd,ladd)
c
      return
c
c*********************************** hii *******************************
c
      entry hii(hiim,isym,numi)
c
      lad=arr+kadd(os(asm+1)+mn(isym+1))
c
      if (abs(a).lt.acrcy) go to 110
      if (abs(b).lt.acrcy) go to 120
c
      jisv=1
      do 102 i=1,numi
         ji=jisv
         ij=i
cdir$ ivdep
         do 101 j=1,i-1
            t=b*int(lad+intb)
            hiim(ij)=t+a*int(lad+intal)
            hiim(ji)=t+a*int(lad+intau)
            lad=lad+3
            ij=ij+numi
            ji=ji+1
  101    continue
         hiim(ij)=b*int(lad+intbd)+a*int(lad+intad)
         lad=lad+3
         jisv=jisv+numi
  102 continue
      return
c
  110 continue
      jisv=1
      do 112 i=1,numi
         ji=jisv
         ij=i
cdir$ ivdep
         do 111 j=1,i-1
            t=b*int(lad+intb)
            hiim(ij)=t
            hiim(ji)=t
            ji=ji+1
            ij=ij+numi
            lad=lad+3
  111    continue
         hiim(ij)=b*int(lad+intbd)
         jisv=jisv+numi
         lad=lad+3
  112 continue
      return
c
  120 continue
      jisv=1
      do 122 i=1,numi
         ji=jisv
         ij=i
cdir$ ivdep
         do 121 j=1,i-1
            hiim(ij)=a*int(lad+intal)
            hiim(ji)=a*int(lad+intau)
            ji=ji+1
            ij=ij+numi
            lad=lad+3
  121    continue
         hiim(ij)=a*int(lad+intad)
         jisv=jisv+numi
         lad=lad+3
  122 continue
c
      return
c
c************************************ hiis *****************************
c
      entry hiis(hiism,isym,numi)
c
      lad=arr+kadd(os(asm+1)+mn(isym+1))
c
      if (abs(a).lt.acrcy) go to 210
      if (abs(b).lt.acrcy) go to 220
c
      jisv=1
      do 202 i=1,numi
         ji=jisv
         ij=i
cdir$ ivdep
         do 201 j=1,i-1
            t=b*int(lad+intb )+a*int(lad+intal)
            hiism(ij)=t
            hiism(ji)=t
            ji=ji+1
            ij=ij+numi
            lad=lad+3
  201    continue
         hiism(ij)=0.0
         jisv=jisv+numi
         lad=lad+3
  202 continue
c
      return
c
  210 continue
      jisv=1
      do 212 i=1,numi
         ji=jisv
         ij=i
cdir$ ivdep
         do 211 j=1,i-1
            t=b*int(lad+intb)
            hiism(ij)=t
            hiism(ji)=t
            ji=ji+1
            ij=ij+numi
            lad=lad+3
  211    continue
         hiism(ij)=0.0
         jisv=jisv+numi
         lad=lad+3
  212 continue
c
      return
c
  220 continue
      jisv=1
cdir$ ivdep
      do 222 i=1,numi
         ji=jisv
         ij=i
         do 221 j=1,i-1
            t=a*int(lad+intal)
            hiism(ij)=t
            hiism(ji)=t
            ji=ji+1
            ij=ij+numi
            lad=lad+3
  221    continue
         hiism(ij)=0.0
         jisv=jisv+numi
         lad=lad+3
  222 continue
c
      return
c
c********************************** hij ********************************
c
      entry hij(hijm,isym,jsym,numi,numj)
c
      lad=arr+kadd(os(asm+1)+mn(isym+1))
c
      if (abs(a).lt.acrcy) go to 310
      if (abs(b).lt.acrcy) go to 320
c
      do 302 i=1,numi
         ij=i
         do 301 j=1,numj
            hijm(ij)=a*int(lad+intal)+b*int(lad+intb)
            ij=ij+numi
            lad=lad+3
  301    continue
  302 continue
c
      return
c
  310 continue
      lad=lad+intb
      do 312 i=1,numi
         ij=i
         do 311 j=1,numj
            hijm(ij)=b*int(lad)
            ij=ij+numi
            lad=lad+3
  311    continue
  312 continue
c
      return
c
  320 continue
      lad=lad+intal
      do 322 i=1,numi
         ij=i
         do 321 j=1,numj
            hijm(ij)=a*int(lad)
            ij=ij+numi
            lad=lad+3
  321    continue
  322 continue
      return
c
c******************************* hi ************************************
c
      entry hi(him,isym,numi)
c
      lad=arr+ladd(os(asm+1)+mn(isym+1))
      do 400 i=1,numi
         him(i)=a*int(lad+intal)+b*int(lad+intb)
         lad=lad+3
  400 continue
      return
c
c******************************* hs ************************************
c
      entry hs(hsm,isym,numi)
c
      a2=a*sqrt(2.0d+00)
      lad=arr+ladd(os(asm+1)+mn(isym+1))+2
      do 410 i=1,numi
         hsm(i)=a2*int(lad)
         lad=lad+3
  410 continue
      return
c
c
      end
