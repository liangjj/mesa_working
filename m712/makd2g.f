*deck @(#)makd2g.f	5.1  11/6/94
      subroutine makd2g(g,dg,d2g,b00,b10,bp01,c00,cp00,dc00,dcp00,
     #                 nmax,mmax,nv,lenv,ncentr,nd2,cen)
c
c
c***begin prologue     makd2g
c***date written       851106   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           derivative integrals, rys polynomials,
c                      two-center integrals
c
c***author             saxe, paul,    (lanl)
c***source             @(#)makd2g.f	5.1   11/6/94
c***purpose            to form the derivatives of the two-center
c                      integrals g used in derivatives.
c
c***description        #
c
c
c***references         (none)
c
c***routines called    (none)
c
c   common blocks      (none)
c
c***end prologue       makd2g
c
      implicit integer (a-z)
c
      real*8 g(nv,3,0:nmax,0:mmax)
      real*8 dg(nv,3,0:nmax,0:mmax,ncentr)
      real*8 d2g(nv,3,0:nmax,0:mmax,nd2)
      real*8 b00(nv)
      real*8 b10(nv)
      real*8 bp01(nv)
      real*8 c00(lenv,3)
      real*8 cp00(lenv,3)
      real*8 dc00(lenv,ncentr)
      real*8 dcp00(lenv,ncentr)
      real*8 one
      real*8 scalr1
      real*8 scalr2
      integer cen(2,nd2)
      integer n
      integer m
      integer coord
      integer c1
      integer c2
      integer nmax
      integer mmax
c
      data one /1.0d+00/
      save one
c
c     ----- g(0,0) -----
c
      do 21 d2=1,nd2
         call rzero(d2g(1,1,0,0,d2),nv*3)
   21 continue
c
      if (nmax.eq.0.and.mmax.eq.0) go to 1000
c
c     ----- g(1,0)=c00 -----
c
      if (nmax.gt.0) then
         do 2 d2=1,nd2
            call rzero(d2g(1,1,1,0,d2),nv*3)
    2    continue
      end if
c
c     ----- g(0,1)=cp00 -----
c
      if (mmax.gt.0) then
         do 4 d2=1,nd2
            call rzero(d2g(1,1,0,1,d2),nv*3)
    4    continue
      end if
c
c     ----- g(1,1)=b00*g(0,0)+cp00*g(1,0) -----
c
      if (nmax.gt.0.and.mmax.gt.0) then
         do 7 d2=1,nd2
            c1=cen(1,d2)
            c2=cen(2,d2)
            do 6 coord=1,3
               do 5 i=1,nv
                  d2g(i,coord,1,1,d2)=
     $                 dcp00(i,c1)*dg(i,coord,1,0,c2)+
     #                 dcp00(i,c2)*dg(i,coord,1,0,c1)
    5          continue
    6       continue
    7    continue
      end if
c
c     ----- g(n,0)=(n-1)b10*g(n-2,0)+c00*g(n-1,0) -----
c
      do 13 n=2,nmax
         scalr1=n-1
         do 12 d2=1,nd2
            c1=cen(1,d2)
            c2=cen(2,d2)
            do 11 coord=1,3
               do 10 i=1,nv
                  d2g(i,coord,n,0,d2)=
     #                 scalr1*b10(i)*d2g(i,coord,n-2,0,d2)+
     #                 dc00(i,c1)*dg(i,coord,n-1,0,c2)+
     $                 dc00(i,c2)*dg(i,coord,n-1,0,c1)+
     #                 c00(i,coord)*d2g(i,coord,n-1,0,d2)
   10          continue
   11       continue
   12    continue
   13 continue
c
c     ----- g(n,1)=g(n-1,0)*b00*n+cp00*g(n,0) -----
c
      if (mmax.gt.0) then
         do 50 n=2,nmax
            scalr1=n
            do 43 d2=1,nd2
               c1=cen(1,d2)
               c2=cen(2,d2)
               do 42 coord=1,3
                  do 41 i=1,nv
                     d2g(i,coord,n,1,d2)=
     #                    scalr1*b00(i)*d2g(i,coord,n-1,0,d2)+
     #                    dcp00(i,c1)*dg(i,coord,n,0,c2)+
     #                    dcp00(i,c2)*dg(i,coord,n,0,c1)+
     #                    cp00(i,coord)*d2g(i,coord,n,0,d2)
   41             continue
   42          continue
   43       continue
   50    continue
      end if
c
c     ----- g(0,m)=g(0,m-2)*bp01*(m-1)+cp00*g(0,m-1) -----
c
      do 70 m=2,mmax
         scalr1=m-1
         do 63 d2=1,nd2
            c1=cen(1,d2)
            c2=cen(2,d2)
            do 62 coord=1,3
               do 61 i=1,nv
                  d2g(i,coord,0,m,d2)=
     #                 scalr1*bp01(i)*d2g(i,coord,0,m-2,d2)+
     #                 dcp00(i,c1)*dg(i,coord,0,m-1,c2)+
     #                 dcp00(i,c2)*dg(i,coord,0,m-1,c1)+
     #                 cp00(i,coord)*d2g(i,coord,0,m-1,d2)
   61          continue
   62       continue
   63    continue
   70 continue
c
c     ----- g(1,m)=g(0,m-1)*b00*m+c00*g(0,m) -----
c
      if (nmax.gt.0) then
         do 90 m=2,mmax
            scalr1=m
            do 83 d2=1,nd2
               c1=cen(1,d2)
               c2=cen(2,d2)
               do 82 coord=1,3
                  do 81 i=1,nv
                     d2g(i,coord,1,m,d2)=
     #                    scalr1*b00(i)*d2g(i,coord,0,m-1,d2)+
     #                    dc00(i,c1)*dg(i,coord,0,m,c2)+
     #                    dc00(i,c2)*dg(i,coord,0,m,c1)+
     #                    c00(i,coord)*d2g(i,coord,0,m,d2)
   81             continue
   82          continue
   83       continue
   90    continue
      end if
c
c     ----- g(n,m)=g(n-2,m)*b10*(n-1)+g(n-1,m-1)*b00*m+c00*g(n-1,m) --
c
      if (nmax.gt.1.and.mmax.gt.1) then
         do 120 m=2,mmax
            scalr1=m
            do 110 n=2,nmax
               scalr2=n-1
               do 103 d2=1,nd2
                  c1=cen(1,d2)
                  c2=cen(2,d2)
                  do 102 coord=1,3
                     do 101 i=1,nv
                        d2g(i,coord,n,m,d2)=
     #                       scalr2*b10(i)*d2g(i,coord,n-2,m,d2)+
     #                       scalr1*b00(i)*d2g(i,coord,n-1,m-1,d2)+
     #                       dc00(i,c1)*dg(i,coord,n-1,m,c2)+
     #                       dc00(i,c2)*dg(i,coord,n-1,m,c1)+
     #                       c00(i,coord)*d2g(i,coord,n-1,m,d2)
  101                continue
  102             continue
  103          continue
  110       continue
  120    continue
      end if
c
 1000 continue
c
c
      return
      end
