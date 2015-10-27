*deck rkck
      subroutine rkck(y,dy,yout,yerr,t,x,h,n)
      implicit integer(a-z)
      real*8 h, x, dy, y, yerr, yout, t, dum
      dimension y(n), dy(n), yout(n), yerr(n), t(n,6)
      real*8 a2, a3, a4, a5, a6
      real*8 b21, b31, b32, b41, b42, b43, b51, b52, b53
      real*8 b54, b61, b62, b63, b64, b65
      real*8 c1, c3, c4, c6
      real*8 dc1 ,dc3, dc4, dc5, dc6
      parameter ( a2=.2d0, a3=.3d0, a4=.6d0, a5=1.d0, a6=.875d0 )
      parameter ( b21=.2d0, b31=.075d0, b32=.225d0, b41=.3d0, 
     1            b42=-.9d0, b43=1.2d0, b51=-.2037037037037037037037d0,
     2            b52=2.5d0, b53=-2.59259259259259259259259d0,
     3            b54=1.296296296296296296296296d0,
     4            b61=.029495804398148148148148148d0,
     5            b62=.341796875d0, b63=.041594328703703703703703d0,
     6            b64=.400345413773148148148148d0,
     7            b65=.061767578125d0)
      parameter( c1=.0978835978835978835978835978835d0,
     1           c3=.402576489533011272141706924315d0,
     2           c4=.2104377104377104377104377104377d0,
     3           c6=.289102202145680406549971767363d0)
      parameter ( dc1=c1 -.102177372685185185185185185185d0,
     1            dc3=c3 -.383907903439153439153439153439d0,
     2            dc4=c4 -.244597160683606112668414865720d0,
     3            dc5=-.019321986607142857142857142857d0,
     4            dc6=c6 -.25d0 )
      do 10 i=1,n
        t(i,6)=y(i)+b21*h*dy(i)
 10   continue
      call df(x+a2*h,t(1,6),t(1,1),dum,idum,n)
      do 20 i=1,n
        t(i,6)=y(i)+h*(b31*dy(i)+b32*t(i,1))
 20   continue
      call df(x+a3*h,t(1,6),t(1,2),dum,idum,n)
      do 30 i=1,n
        t(i,6)=y(i)+h*(b41*dy(i)+b42*t(i,1)+b43*t(i,2))
 30   continue
      call df(x+a4*h,t(1,6),t(1,3),dum,idum,n)
      do 40 i=1,n
        t(i,6)=y(i)+h*(b51*dy(i)+b52*t(i,1)+b53*t(i,2)+b54*t(i,3))
 40   continue
      call df(x+a5*h,t(1,6),t(1,4),dum,idum,n)
      do 50 i=1,n
        t(i,6)=y(i)+h*(b61*dy(i)+b62*t(i,1)+b63*t(i,2)+b64*t(i,3)+
     1                   b65*t(i,4))
 50   continue
      call df(x+a6*h,t(1,6),t(1,5),dum,idum,n)
      do 60 i=1,n
        yout(i)=y(i)+h*(c1*dy(i)+c3*t(i,2)+c4*t(i,3)+c6*t(i,5))
 60   continue
      do 70 i=1,n
        yerr(i)=h*(dc1*dy(i)+dc3*t(i,2)+dc4*t(i,3)+dc5*t(i,4)+dc6*
     1             t(i,5))
  70  continue
      return
      end


