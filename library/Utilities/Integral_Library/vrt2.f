*deck @(#)vrt2.f	5.1  11/6/94
      subroutine vrt2(x,u,w,n,switch)
c
c***purpose: to calculate two roots and weights for a rys polynomial
c
      implicit integer (a-z)
c
      character*4 itoc
      real*8 pie4,r12,r22,w22,cf1,cf2,x(n),u(n,2),w(n,2)
c
      common /cfrt/   pie4,r12,r22,w22,cf1(52),cf2(138)
c
      go to (10,20,30,40,50,60,70,80,90), switch
      call lnkerr('unrecognized switch in vrt2:'//itoc(switch))
c
c
   10 continue
c
c     ----- x = 0.0 to 3.0d-07 -----
c
         do 1 i=1,n
            u(i,1)=cf2(1)+cf2(2)*x(i)
            u(i,2)=cf2(3)+cf2(4)*x(i)
            w(i,1)=cf2(5)+cf2(6)*x(i)
            w(i,2)=cf2(7)+cf2(8)*x(i)
    1    continue
         return
c
   20 continue
c
c     ----- x = 3.0d-07 to 1.0d+00 -----
c
         call vpoly(w(1,2),x,cf1(5),10,n)
         call vpoly(u(1,1),x,cf2(9),9,n)
         call vpoly(u(1,2),x,cf2(18),9,n)
         do 21 i=1,n
            w(i,1)=2.0d+00*x(i)*w(i,2)+exp(-x(i))
   21    continue
         go to 100
c
   30 continue
c
c     ----- x = 1.0 to 3.0 -----
c
         do 31 i=1,n
            w(i,1)=x(i)-2.0d+00
   31    continue
         call vpoly(w(1,2),w(1,1),cf1(15),12,n)
         call vpoly(u(1,1),w(1,1),cf2(27),11,n)
         call vpoly(u(1,2),w(1,1),cf2(38),11,n)
         do 32 i=1,n
            w(i,1)=2.0d+00*x(i)*w(i,2)+exp(-x(i))
   32    continue
         go to 100
c
   40 continue
c
c     ----- x = 3.0 to 5.0 -----
c
         do 41 i=1,n
            w(i,1)=x(i)-4.0d+00
   41    continue
         call vpoly(w(1,2),w(1,1),cf1(27),12,n)
         call vpoly(u(1,1),w(1,1),cf2(49),10,n)
         call vpoly(u(1,2),w(1,1),cf2(59),11,n)
         do 42 i=1,n
            w(i,1)=2.0d+00*x(i)*w(i,2)+exp(-x(i))
   42    continue
         go to 100
c
   50 continue
c
c     ----- x = 5.0 to 10.0 -----
c
         do 53 i=1,n
            w(i,2)=x(i)-7.5d+00
   53    continue
         call vpoly(u(1,1),w(1,2),cf2(70),15,n)
         call vpoly(u(1,2),w(1,2),cf2(85),14,n)
         do 51 i=1,n
            w(i,2)=exp(-x(i))
   51    continue
         call vpolyd(w(1,1),x,cf1(39),7,n)
         do 52 i=1,n
            w(i,1)=w(i,1)*w(i,2)+sqrt(pie4/x(i))
            w(i,2)=(w(i,1)-w(i,2))/(2.0d+00*x(i))
   52    continue
         go to 100
c
   60 continue
c
c     ----- x = 10.0 to 15.0 -----
c
         do 61 i=1,n
            w(i,2)=exp(-x(i))
   61    continue
         call vpolyd(w(1,1),x,cf2(104),5,n)
         call vpoly(u(1,1),x,cf2(99),5,n)
         do 62 i=1,n
            u(i,1)=(u(i,1)+w(i,1))*w(i,2)+r12/(x(i)-r12)
   62    continue
         call vpolyd(w(1,1),x,cf2(113),5,n)
         call vpoly(u(1,2),x,cf2(109),4,n)
         do 63 i=1,n
            u(i,2)=(u(i,2)+w(i,1))*w(i,2)+r22/(x(i)-r22)
   63    continue
         call vpolyd(w(1,1),x,cf1(46),4,n)
         do 64 i=1,n
            w(i,1)=w(i,1)*w(i,2)+sqrt(pie4/x(i))
            w(i,2)=(w(i,1)-w(i,2))/(2.0d+00*x(i))
   64    continue
         go to 100
c
   70 continue
c
c     ----- x = 15.0 to 33.0 -----
c
         do 71 i=1,n
            w(i,2)=exp(-x(i))
   71    continue
         call vpolyd(w(1,1),x,cf2(123),3,n)
         call vpoly(u(1,1),x,cf2(118),5,n)
         do 72 i=1,n
            u(i,1)=(u(i,1)+w(i,1))*w(i,2)+r12/(x(i)-r12)
   72    continue
         call vpolyd(w(1,1),x,cf2(130),3,n)
         call vpoly(u(1,2),x,cf2(126),4,n)
         do 73 i=1,n
            u(i,2)=(u(i,2)+w(i,1))*w(i,2)+r22/(x(i)-r22)
   73    continue
         call vpolyd(w(1,1),x,cf1(50),3,n)
         do 74 i=1,n
            w(i,1)=w(i,1)*w(i,2)+sqrt(pie4/x(i))
            w(i,2)=(w(i,1)-w(i,2))/(2.0d+00*x(i))
   74    continue
         go to 100
c
   80 continue
c
c     ----- x = 33.0 to 40.0 -----
c
         do 81 i=1,n
            w(i,1)=exp(-x(i))
   81    continue
         do 82 i=1,n
            u(i,1)=(cf2(133)*x(i)+cf2(134))*w(i,1)+r12/(x(i)-r12)
            u(i,2)=(cf2(135)*x(i)+cf2(136))*w(i,1)+r22/(x(i)-r22)
            w(i,2)=(cf2(137)*x(i)+cf2(138))*w(i,1)
   82    continue
         do 83 i=1,n
            w(i,1)=sqrt(pie4/x(i))
            w(i,2)=w(i,2)+w22*w(i,1)
            w(i,1)=w(i,1)-w(i,2)
   83    continue
         return
c
   90 continue
c
c     ----- x = 40.0 to infinity -----
c
         do 91 i=1,n
            w(i,1)=sqrt(pie4/x(i))
            u(i,1)=r12/(x(i)-r12)
            u(i,2)=r22/(x(i)-r22)
            w(i,2)=w22*w(i,1)
            w(i,1)=w(i,1)-w(i,2)
   91    continue
         return
c
c
  100 continue
c
c     ----- the finishing touches -----
c
      do 101 i=1,n
         w(i,2)=((w(i,2)-w(i,1))*u(i,1)+w(i,2))*(1.0d+00+u(i,2))/
     #           (u(i,2)-u(i,1))
         w(i,1)=w(i,1)-w(i,2)
  101 continue
c
c
      return
      end
