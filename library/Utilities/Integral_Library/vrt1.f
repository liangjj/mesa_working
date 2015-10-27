*deck @(#)vrt1.f	5.1  11/6/94
      subroutine vrt1(x,u,w,n,switch,t1,t2)
c
c***module to calculate the roots and weight for a rys polynomial of
c   degree one.
c
      implicit integer (a-z)
c
      character*4 itoc
      real*8 pie4,r12,r22,w22,cf1,cf2,x(n),u(n),w(n),t1(n),t2(n)
c
      common /cfrt/   pie4,r12,r22,w22,cf1(52),cf2(138)
c
      go to (1,2,3,4,5,6,7,8), switch
      call lnkerr ('unrecognized switch in vrt1:'//itoc(switch))
c
c
    1 continue
c
c     ----- x nearly 0.0 -----
c
         do 11 i=1,n
            u(i)=0.5d+00-x(i)/5.0d+00
            w(i)=1.0d+00-x(i)/3.0d+00
   11    continue
         return
c
    2 continue
c
c     ----- x = 0.0 to 1.0 -----
c
         call vpoly(t2,x,cf1(5),10,n)
         go to 100
c
    3 continue
c
c     ----- x = 1.0 to 3.0 -----
c
         do 31 i=1,n
            t1(i)=x(i)-2.0d+00
   31    continue
         call vpoly(t2,t1,cf1(15),12,n)
         go to 100
    4 continue
c
c     ----- x = 3.0 to 5.0 -----
c
         do 41 i=1,n
            t1(i)=x(i)-4.0d+00
   41    continue
         call vpoly(t2,t1,cf1(27),12,n)
         go to 100
    5 continue
c
c     ----- x = 5.0 to 10.0 -----
c
         do 51 i=1,n
            t1(i)=exp(-x(i))
   51    continue
         call vpolyd(t2,x,cf1(39),7,n)
         go to 110
c
    6 continue
c
c     ----- x = 10.0 to 15.0 -----
c
         do 61 i=1,n
            t1(i)=exp(-x(i))
   61    continue
         call vpolyd(t2,x,cf1(46),4,n)
         go to 110
    7 continue
c
c     ----- x = 15.0 to 33.0 -----
c
         do 71 i=1,n
            t1(i)=exp(-x(i))
   71    continue
         call vpolyd(t2,x,cf1(50),3,n)
         go to 110
    8 continue
c
c     ----- x = 33.0 to infinity -----
c
         do 81 i=1,n
            w(i)=sqrt(pie4/x(i))
            u(i)=0.5d+00/(x(i)-0.5d+00)
   81    continue
         return
c
c
  100 continue
      do 101 i=1,n
         w(i)=2.0d+00*x(i)*t2(i)+exp(-x(i))
         u(i)=t2(i)/(w(i)-t2(i))
  101 continue
      return
c
c
  110 continue
      do 111 i=1,n
         w(i)=t2(i)*t1(i)+sqrt(pie4/x(i))
  111 continue
      do 112 i=1,n
         t2(i)=(w(i)-t1(i))/(2.0d+00*x(i))
         u(i)=t2(i)/(w(i)-t2(i))
  112 continue
      return
c
c
      end
