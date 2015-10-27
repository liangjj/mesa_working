*deck @(#)makdg.f	5.1  11/6/94
      subroutine makdg(g,dg,f00,b00,b10,bp01,c00,cp00,dc00,dcp00,
     #                 nmax,mmax,nv,lenv,ncentr)
c
c
c***begin prologue     makdg
c***date written       851106   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           derivative integrals, rys polynomials,
c                      two-center integrals
c
c***author             saxe, paul,    (lanl)
c***source             @(#)makdg.f	5.1   11/6/94
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
c***end prologue       makdg
c
      implicit integer (a-z)
c
      real*8 g(nv,3,0:nmax,0:mmax),dg(nv,3,0:nmax,0:mmax,ncentr)
      real*8 f00(nv),b00(nv),b10(nv),bp01(nv),c00(lenv,3),cp00(lenv,3)
      real*8 dc00(lenv,ncentr),dcp00(lenv,ncentr)
      real*8 one,scalr1,scalr2
      integer n,m,coord
      integer nmax,mmax
c
      data one /1.0d+00/
      save one
c
c     ----- g(0,0) -----
c
      do 21 centr=1,ncentr
         call rzero(dg(1,1,0,0,centr),nv*3)
   21 continue
c
      if (nmax.eq.0.and.mmax.eq.0) go to 1000
c
c     ----- g(1,0)=c00 -----
c
      if (nmax.gt.0) then
         do 2 centr=1,ncentr
                do 1 i=1,nv
                   dg(i,1,1,0,centr)=dc00(i,centr)
                   dg(i,2,1,0,centr)=dc00(i,centr)
                   dg(i,3,1,0,centr)=dc00(i,centr)*f00(i)
    1           continue
    2    continue
      end if
c
c     ----- g(0,1)=cp00 -----
c
      if (mmax.gt.0) then
         do 4 centr=1,ncentr
            do 3 i=1,nv
               dg(i,1,0,1,centr)=dcp00(i,centr)
               dg(i,2,0,1,centr)=dcp00(i,centr)
               dg(i,3,0,1,centr)=dcp00(i,centr)*f00(i)
    3       continue
    4    continue
      end if
c
c     ----- g(1,1)=b00*g(0,0)+cp00*g(1,0) -----
c
      if (nmax.gt.0.and.mmax.gt.0) then
         do 7 centr=1,ncentr
            do 6 coord=1,3
               do 5 i=1,nv
                  dg(i,coord,1,1,centr)=dcp00(i,centr)*g(i,coord,1,0)+
     #                              cp00(i,coord)*dg(i,coord,1,0,centr)
    5          continue
    6       continue
    7    continue
      end if
c
c     ----- g(n,0)=(n-1)b10*g(n-2,0)+c00*g(n-1,0) -----
c
      do 13 n=2,nmax
         scalr1=n-1
         do 12 centr=1,ncentr
            do 11 coord=1,3
               do 10 i=1,nv
                  dg(i,coord,n,0,centr)=
     #                          scalr1*b10(i)*dg(i,coord,n-2,0,centr)+
     #                                 dc00(i,centr)*g(i,coord,n-1,0)+
     #                           c00(i,coord)*dg(i,coord,n-1,0,centr)
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
            do 43 centr=1,ncentr
               do 42 coord=1,3
                  do 41 i=1,nv
                     dg(i,coord,n,1,centr)=
     #                         scalr1*b00(i)*dg(i,coord,n-1,0,centr)+
     #                         dcp00(i,centr)*g(i,coord,n,0)+
     #                         cp00(i,coord)*dg(i,coord,n,0,centr)
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
         do 63 centr=1,ncentr
            do 62 coord=1,3
               do 61 i=1,nv
                  dg(i,coord,0,m,centr)=
     #                        scalr1*bp01(i)*dg(i,coord,0,m-2,centr)+
     #                        dcp00(i,centr)*g(i,coord,0,m-1)+
     #                        cp00(i,coord)*dg(i,coord,0,m-1,centr)
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
            do 83 centr=1,ncentr
               do 82 coord=1,3
                  do 81 i=1,nv
                     dg(i,coord,1,m,centr)=
     #                           scalr1*b00(i)*dg(i,coord,0,m-1,centr)+
     #                           dc00(i,centr)*g(i,coord,0,m)+
     #                           c00(i,coord)*dg(i,coord,0,m,centr)
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
               do 103 centr=1,ncentr
                  do 102 coord=1,3
                     do 101 i=1,nv
                        dg(i,coord,n,m,centr)=
     #                      scalr2*b10(i)*dg(i,coord,n-2,m,centr)+
     #                      scalr1*b00(i)*dg(i,coord,n-1,m-1,centr)+
     #                      dc00(i,centr)*g(i,coord,n-1,m)+
     #                      c00(i,coord)*dg(i,coord,n-1,m,centr)
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
