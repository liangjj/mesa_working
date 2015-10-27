*deck @(#)prims.f	5.1  11/6/94
      subroutine prims(a,ar,xyza,xyzam1,xyzam3,nij,
     #                 b,br,xyzb,xyzbm1,xyzbm3,nkl,minkl,maxkl,
     #                 f00,b00,b10,bp01,c00,cp00,
     #                 ab,aplusb,urho,wt,denom,test,
     #                 t1,t2,t3,t4,t5,t6,t7,t8,t9,
     #                 nijkl,lenv,nmax,mmax,nroots,nbatch,
     #                 cutexp,pi252)
c
c
      implicit integer (a-z)
c
      real*8 a(nij),ar(nij),xyza(nij,3),xyzam1(nij,3),xyzam3(nij,3)
      real*8 b(nkl),br(nkl),xyzb(nkl,3),xyzbm1(nkl,3),xyzbm3(nkl,3)
      real*8 f00(lenv),b00(lenv),b10(lenv),bp01(lenv)
      real*8 c00(lenv,3),cp00(lenv,3)
      real*8 ab(nijkl),aplusb(nijkl),urho(lenv),wt(lenv),denom(lenv)
      real*8 t1(nijkl),t2(nijkl),t3(nijkl),t4(lenv),t5(lenv)
      real*8 t6(nijkl),t7(nijkl),t8(nijkl),t9(nijkl)
      real*8 cutexp,pi252
      real*8 scal1,scal2,scal3
      integer test(nijkl)
      integer inp,iout
c
      common /io/     inp,iout
c
c     ------ start timing -----
c
c
c     ----- form the exponential prefactor -----
c
      ijkl=0
      do 2 kl=minkl,maxkl
         scal1=br(kl)
         do 1 ij=1,nij
            ijkl=ijkl+1
            t1(ijkl)=ar(ij)+scal1
    1    continue
    2 continue
c
c     ----- form a test vector of exponent minus cutoff -----
c
      ipt=0
      do 3 ijkl=1,nijkl
         if (t1(ijkl).lt.cutexp) then
            ipt=ipt+1
            test(ipt)=ijkl
         end if
c         test(ijkl)=t1(ijkl)-cutexp
    3 continue
      nbatch=ipt
c
c     ----- and extract the prefactors large enough -----
c
c      call kmprsm(nijkl,test,1,t1,1,t2,1,nbatch)
      do 201 ipt=1,nbatch
         t2(ipt)=t1(test(ipt))
  201 continue
c
c     ----- form ab and a+b -----
c
      ijkl=0
      do 5 kl=minkl,maxkl
         scal1=b(kl)
         do 4 ij=1,nij
            ijkl=ijkl+1
            t1(ijkl)=a(ij)*scal1
            t3(ijkl)=a(ij)+scal1
    4    continue
    5 continue
c
c      call kmprsm(nijkl,test,1,t1,1,ab,1,junk)
c      call kmprsm(nijkl,test,1,t3,1,aplusb,1,junk)
      do 102 ipt=1,nbatch
         ab(ipt)=t1(test(ipt))
         aplusb(ipt)=t3(test(ipt))
  102 continue
c
c     ----- form rho=ab/(a+b) -----
c
      do 6 i=1,nbatch
         t9(i)=ab(i)/aplusb(i)
    6 continue
c
c     ----- form the argument for the rys polynomials -----
c
      ijkl=0
      do 8 kl=minkl,maxkl
         scal1=xyzb(kl,1)
         scal2=xyzb(kl,2)
         scal3=xyzb(kl,3)
         do 7 ij=1,nij
            ijkl=ijkl+1
            t4(ijkl)=(xyza(ij,1)-scal1)**2+(xyza(ij,2)-scal2)**2+
     #               (xyza(ij,3)-scal3)**2
    7    continue
    8 continue
c
c      call kmprsm(nijkl,test,1,t4,1,t3,1,junk)
      do 103 ipt=1,nbatch
         t3(ipt)=t4(test(ipt))
  103 continue
c
      do 9 i=1,nbatch
         t1(i)=t9(i)*t3(i)
    9 continue
c
c     ----- and find the roots -----
c
      call roots(nroots,t1,urho,wt,nbatch,t3,t4,t5,t6,t7,t8)
c
c     ----- and the exponential prefactor -----
c
      do 11 i=1,nbatch
         t2(i)=exp(-t2(i))*pi252/(ab(i)*sqrt(aplusb(i)))
   11 continue
c
      do 13 root=1,nroots
         do 12 i=1,nbatch
            f00(i*nroots-nroots+root)=t2(i)*wt(i*nroots-nroots+root)
   12    continue
   13 continue
c
c     ----- if [ss;ss] block, f00 is the integrals -----
c
      if (nmax+mmax.eq.0) go to 1000
c
c     ----- form u*rho -----
c
      do 101 root=1,nroots
         do 10 i=1,nbatch
            urho(i*nroots-nroots+root)=urho(i*nroots-nroots+root)*t9(i)
   10    continue
  101 continue
c
c     ----- form the denominator for subsequent steps -----
c
      do 15 root=1,nroots
         do 14 i=1,nbatch
            denom(i*nroots-nroots+root)=1.0d+00/(ab(i)+aplusb(i)*
     #                            urho(i*nroots-nroots+root))
   14    continue
   15 continue
c
c
      do 30 coord=1,3
c
c     ----- form c00 -----
c
         if (nmax.gt.0) then
            ijkl=0
            do 23 kl=minkl,maxkl
               scal1=xyzbm1(kl,coord)
               scal2=b(kl)
               do 22 ij=1,nij
                  ijkl=ijkl+1
                  t3(ijkl)=xyzam1(ij,coord)+scal1
                  t4(ijkl)=xyzam1(ij,coord)*scal2
   22          continue
   23       continue
c
c            call kmprsm(nijkl,test,1,t3,1,t1,1,junk)
c            call kmprsm(nijkl,test,1,t4,1,t2,1,junk)
            do 105 ipt=1,nbatch
               t1(ipt)=t3(test(ipt))
               t2(ipt)=t4(test(ipt))
  105       continue
c
            do 25 root=1,nroots
               do 24 i=1,nbatch
                  c00(i*nroots-nroots+root,coord)=(t2(i)+
     #                         t1(i)*urho(i*nroots-nroots+root))*
     #                         denom(i*nroots-nroots+root)
   24          continue
   25       continue
         end if
c
c     ----- form cp00 -----
c
         if (mmax.gt.0) then
            ijkl=0
            do 27 kl=minkl,maxkl
               scal1=xyzbm3(kl,coord)
               do 26 ij=1,nij
                  ijkl=ijkl+1
                  t3(ijkl)=xyzam3(ij,coord)+scal1
                  t4(ijkl)=a(ij)*scal1
   26          continue
   27       continue
c
c            call kmprsm(nijkl,test,1,t3,1,t1,1,junk)
c            call kmprsm(nijkl,test,1,t4,1,t2,1,junk)
            do 106 ipt=1,nbatch
               t1(ipt)=t3(test(ipt))
               t2(ipt)=t4(test(ipt))
  106       continue
c
            do 29 root=1,nroots
               do 28 i=1,nbatch
                  cp00(i*nroots-nroots+root,coord)=(t2(i)+
     #                         t1(i)*urho(i*nroots-nroots+root))*
     #                         denom(i*nroots-nroots+root)
   28          continue
   29       continue
         end if
   30 continue
c
c     ----- double the denominator for b00, b10 and bp01 -----
c
      do 31 i=1,nbatch*nroots
         denom(i)=denom(i)*0.5d+00
   31 continue
c
c     ----- form b00 -----
c
      if (nmax.gt.0.and.mmax.gt.0) then
         do 33 i=1,nbatch*nroots
            b00(i)=urho(i)*denom(i)
   33    continue
      end if
c
c     ----- form b10 -----
c
      if (nmax.gt.1) then
         ijkl=0
         do 51 kl=minkl,maxkl
            scal1=b(kl)
            do 50 ij=1,nij
               ijkl=ijkl+1
               t1(ijkl)=scal1
   50       continue
   51    continue
c
c         call kmprsm(nijkl,test,1,t1,1,t2,1,junk)
         do 107 ipt=1,nbatch
            t2(ipt)=t1(test(ipt))
  107    continue
c
         do 36 root=1,nroots
            do 35 i=1,nbatch
               b10(i*nroots-nroots+root)=(t2(i)+
     #                         urho(i*nroots-nroots+root))*
     #                         denom(i*nroots-nroots+root)
   35       continue
   36    continue
      end if
c
c     ----- form bp01 -----
c
      if (mmax.gt.1) then
         ijkl=0
         do 53 kl=minkl,maxkl
            do 52 ij=1,nij
               ijkl=ijkl+1
               t1(ijkl)=a(ij)
   52       continue
   53    continue
c
c         call kmprsm(nijkl,test,1,t1,1,t2,1,junk)
         do 108 ipt=1,nbatch
            t2(ipt)=t1(test(ipt))
  108    continue
c
         do 39 root=1,nroots
            do 38 i=1,nbatch
               bp01(i*nroots-nroots+root)=(t2(i)+
     #                         urho(i*nroots-nroots+root))*
     #                         denom(i*nroots-nroots+root)
   38       continue
   39    continue
      end if
c
c     ----- stop timing -----
c
 1000 continue
c
c
      return
      end
