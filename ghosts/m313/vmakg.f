*deck @(#)vmakg.f	1.1  11/30/90
      subroutine vmakg(g,f00,b00,b10,bp01,c00,cp00,nmax,mmax,nv,lenv)
c
c     module to form the two-centre two-dimensional integrals g.
c                              paul saxe  23 april 1984  brl
c
c     vector version.
c  paul saxe                 12 july 1984              lanl
c
      implicit integer (a-z)
c
      real*8 g(nv,0:nmax,0:mmax,3)
      real*8 f00(nv),b00(nv),b10(nv),bp01(nv),c00(lenv,3),cp00(lenv,3)
      real*8 one,scalr1,scalr2
      integer n,m,coord
      integer nmax,mmax
c
      data one /1.0d+00/
c
c     ----- g(0,0) -----
c
      call vfill(g(1,0,0,1),one,nv)
      call vfill(g(1,0,0,2),one,nv)
      call vmove(g(1,0,0,3),f00,nv)
c
      if (nmax.eq.0.and.mmax.eq.0) go to 1000
c
c     ----- g(1,0)=c00 -----
c
      if (nmax.gt.0) then
         call vmove(g(1,1,0,1),c00(1,1),nv)
         call vmove(g(1,1,0,2),c00(1,2),nv)
         call vmul(g(1,1,0,3),c00(1,3),f00,nv)
      end if
c
c     ----- g(0,1)=cp00 -----
c
      if (mmax.gt.0) then
         call vmove(g(1,0,1,1),cp00(1,1),nv)
         call vmove(g(1,0,1,2),cp00(1,2),nv)
         call vmul(g(1,0,1,3),cp00(1,3),f00,nv)
      end if
c
c     ----- g(1,1)=b00*g(0,0)+cp00*g(1,0) -----
c
      if (nmax.gt.0.and.mmax.gt.0) then
         do 10 coord=1,3
            do 9 i=1,nv
               g(i,1,1,coord)=b00(i)*g(i,0,0,coord)+
     #                        cp00(i,coord)*g(i,1,0,coord)
    9       continue
   10    continue
      end if
c
c     ----- g(n,0)=(n-1)b10*g(n-2,0)+c00*g(n-1,0) -----
c
      do 30 n=2,nmax
         scalr1=n-1
c        call smul(t1,b10,scalr1,nv)
         do 20 coord=1,3
            do 19 i=1,nv
               g(i,n,0,coord)=scalr1*b10(i)*g(i,n-2,0,coord)+
     #                        c00(i,coord)*g(i,n-1,0,coord)
   19       continue
   20    continue
   30 continue
c
c     ----- g(n,1)=g(n-1,0)*b00*n+cp00*g(n,0) -----
c
      if (mmax.gt.0) then
         do 50 n=2,nmax
            scalr1=n
c           call smul(t1,b00,scalr1,nv)
            do 40 coord=1,3
               do 39 i=1,nv
                  g(i,n,1,coord)=scalr1*b00(i)*g(i,n-1,0,coord)+
     #                           cp00(i,coord)*g(i,n,0,coord)
   39          continue
   40       continue
   50    continue
      end if
c
c     ----- g(0,m)=g(0,m-2)*bp01*(m-1)+cp00*g(0,m-1) -----
c
      do 70 m=2,mmax
         scalr1=m-1
c        call smul(t1,bp01,scalr1,nv)
         do 60 coord=1,3
            do 59 i=1,nv
               g(i,0,m,coord)=scalr1*bp01(i)*g(i,0,m-2,coord)+
     #                        cp00(i,coord)*g(i,0,m-1,coord)
   59       continue
   60    continue
   70 continue
c
c     ----- g(1,m)=g(0,m-1)*b00*m+c00*g(0,m) -----
c
      if (nmax.gt.0) then
         do 90 m=2,mmax
            scalr1=m
c           call smul(t1,b00,scalr1,nv)
            do 80 coord=1,3
               do 79 i=1,nv
                  g(i,1,m,coord)=scalr1*b00(i)*g(i,0,m-1,coord)+
     #                           c00(i,coord)*g(i,0,m,coord)
   79          continue
   80       continue
   90    continue
      end if
c
c     ----- g(n,m)=g(n-2,m)*b10*(n-1)+g(n-1,m-1)*b00*m+c00*g(n-1,m) --
c
      if (nmax.gt.1.and.mmax.gt.1) then
         do 120 m=2,mmax
            scalr1=m
c           call smul(t1,b00,scalr1,nv)
            do 110 n=2,nmax
               scalr2=n-1
c              call smul(t2,b10,scalr2,nv)
               do 100 coord=1,3
                  do 99 i=1,nv
                     g(i,n,m,coord)=scalr2*b10(i)*g(i,n-2,m,coord)+
     #                              scalr1*b00(i)*g(i,n-1,m-1,coord)+
     #                              c00(i,coord)*g(i,n-1,m,coord)
   99             continue
  100          continue
  110       continue
  120    continue
      end if
c
 1000 continue
c
c
      return
      end
