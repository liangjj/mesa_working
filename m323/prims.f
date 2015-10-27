*deck @(#)prims.f	5.1  11/6/94
      subroutine prims(aij,xyza,xyzam1,xyzam3,nij,
     #                 bkl,xyzb,xyzbm1,xyzbm3,nkl,
     #                 f00,b00,b10,bp01,c00,cp00,
     #                 ab,aplusb,urho,wt,denom,a,b,rho,expon,
     #                 t1,t2,t3,t4,t5,t6,t7,t8,
     #                 nv,lenv,nmax,mmax,nroots,
     #                 index,len)
c
c 17 november 1985 extensively modified by pws at lanl
c    original version from m312 practically rewritten for m712
c
      implicit integer (a-z)
c
      real*8 aij(nij),xyza(nij,3),xyzam1(nij,3),xyzam3(nij,3)
      real*8 bkl(nkl),xyzb(nkl,3),xyzbm1(nkl,3),xyzbm3(nkl,3)
      real*8 f00(nroots,nv),b00(nroots,nv),b10(nroots,nv)
      real*8 bp01(nroots,nv)
      real*8 c00(nroots,nv,3),cp00(nroots,nv,3)
      real*8 ab(nv),aplusb(nv),urho(nroots,nv),wt(nroots,nv)
      real*8 denom(nroots,nv),a(nv),b(nv),rho(nv),expon(nv)
      real*8 t1(nv),t2(nv),t3(nv),t4(lenv),t5(lenv)
      real*8 t6(nv),t7(nv),t8(nv)
      integer index(len,6)
c
c     ----- form ab and a+b -----
c
      do 1 i=1,nv
         a(i)=aij(index(i,5))
         b(i)=bkl(index(i,6))
    1 continue
      do 2 i=1,nv
         ab(i)=a(i)*b(i)
         aplusb(i)=a(i)+b(i)
    2 continue
c
c     ----- form rho=ab/(a+b) -----
c
      do 6 i=1,nv
         rho(i)=ab(i)/aplusb(i)
    6 continue
c
c     ----- form the argument for the rys polynomials -----
c
      do 7 i=1,nv
         t2(i)=(xyza(index(i,5),1)-xyzb(index(i,6),1))**2+
     #         (xyza(index(i,5),2)-xyzb(index(i,6),2))**2+
     #         (xyza(index(i,5),3)-xyzb(index(i,6),3))**2
    7 continue
c
      do 9 i=1,nv
         t1(i)=rho(i)*t2(i)
    9 continue
c
c     ----- and find the roots -----
c
      call roots(nroots,t1,urho,wt,nv,t3,t4,t5,t6,t7,t8)
c
      do 13 root=1,nroots
         do 12 i=1,nv
            f00(root,i)=expon(i)*wt(root,i)
   12    continue
   13 continue
c
c     ----- form u*rho -----
c
      do 101 root=1,nroots
         do 10 i=1,nv
            urho(root,i)=urho(root,i)*rho(i)
   10    continue
  101 continue
c
c     ----- if [ss;ss] block, f00 is the integrals -----
c
      if (nmax+mmax.eq.0) go to 1000
c
c     ----- form the denominator for subsequent steps -----
c
      do 15 root=1,nroots
         do 14 i=1,nv
            denom(root,i)=1.0d+00/(ab(i)+aplusb(i)*
     #                            urho(root,i))
   14    continue
   15 continue
c
c
      do 30 coord=1,3
c
c     ----- form c00 -----
c
         if (nmax.gt.0) then
            do 25 root=1,nroots
               do 24 i=1,nv
                  c00(root,i,coord)=(xyzam1(index(i,5),coord)*b(i)+
     #            (xyzam1(index(i,5),coord)+xyzbm1(index(i,6),coord))*
     #                          urho(root,i))*denom(root,i)
   24          continue
   25       continue
         end if
c
c     ----- form cp00 -----
c
         if (mmax.gt.0) then
            do 29 root=1,nroots
               do 28 i=1,nv
                  cp00(root,i,coord)=(xyzbm3(index(i,6),coord)*a(i)+
     #             (xyzam3(index(i,5),coord)+xyzbm3(index(i,6),coord))*
     #                               urho(root,i))*denom(root,i)
   28          continue
   29       continue
         end if
   30 continue
c
c     ----- double the denominator for b00, b10 and bp01 -----
c
      do 70 root=1,nroots
         do 31 i=1,nv
            denom(root,i)=denom(root,i)*0.5d+00
   31    continue
   70 continue
c
c     ----- form b00 -----
c
      if (nmax.gt.0.and.mmax.gt.0) then
         do 71 root=1,nroots
            do 33 i=1,nv
               b00(root,i)=urho(root,i)*denom(root,i)
   33       continue
   71    continue
      end if
c
c     ----- form b10 -----
c
      if (nmax.gt.1) then
         do 36 root=1,nroots
            do 35 i=1,nv
               b10(root,i)=(b(i)+urho(root,i))*denom(root,i)
   35       continue
   36    continue
      end if
c
c     ----- form bp01 -----
c
      if (mmax.gt.1) then
         do 39 root=1,nroots
            do 38 i=1,nv
               bp01(root,i)=(a(i)+urho(root,i))*denom(root,i)
   38       continue
   39    continue
      end if
c
c
 1000 continue
c
c
c
      return
      end
