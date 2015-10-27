*deck bffgh
      subroutine bffgh (z,l,ba,bjp,bjpp,bb,byp,bypp,qsq,ikq,y)
c     c3 ucsd bffgh
      implicit real*8 (a-h,o-z)
      dimension bj(5), by(5), r(3)
c
c     revised 1-22-68
c
      if (l.gt.0) go to 20
      if (ikq.ge.2.and.qsq.lt.0.d0) go to 10
      ba=sin(z)/z
      bb=-cos(z)/z
      uta=ba/z
      utb=bb/z
      bjp=-bb-uta
      byp=ba-utb
      ujtemp=(2.0d0/z**2)-1.0d0
      bjpp=ba*ujtemp+2.0d0*utb
      bypp=bb*ujtemp-2.0d0*uta
      return
   10 exx=exp(z)
      exxi=1.d0/exx
      ba=(exx-exxi)*0.5d0
      ba=ba/z
      bb=exxi/z
      bbb=(exx+exxi)*0.5d0
      bbb=bbb/z
      bjp=bbb-ba/z
      byp=-bb-bb/z
      bjpp=ba-2.d0*bjp/z
      bypp=bb-2.d0*byp/z
      return
   20 ysave=y
      if (y-1.d0) 40,40,30
   30 y=1.
   40 fl=l
      l1=l+1
      zb=z
      ze=0.d0
      fl1=fl+1.d0
      zsq=z*z
      flr=sqrt(fl*fl1)
      sine=sin(z)
      cose=cos(z)
      if (ikq-1) 60,60,70
   50 if (z-flr) 410,410,570
   60 if (z-1.d0) 200,90,90
   70 if (qsq) 80,80,60
   80 if (z-1.d0) 330,320,320
   90 bj(1)=sine/zb
      sign=1.d0
      er=-y
      by(1)=er*cose
  100 by(1)=by(1)/zb
      er=bj(1)/zb
      if (sign) 110,120,120
  110 er=-er
  120 bj(2)=(by(1)/y+er)/y
      by(2)=(by(1)/z-y*bj(1))*y
      flx=2.d0
  130 temp=(2.d0*flx-1.d0)/z
      bj(3)=(temp*bj(2)-bj(1)/y)/y
      if (sign) 140,150,150
  140 bj(3)=-bj(3)
  150 er=temp*by(2)
      by(3)=y*by(1)
      if (sign) 160,170,170
  160 by(3)=-by(3)
  170 by(3)=((er-by(3)))*y
      if (flx-fl1) 190,180,180
  180 if (sign) 340,340,50
  190 bj(1)=bj(2)
      by(1)=by(2)
      bj(2)=bj(3)
      by(2)=by(3)
      flx=flx+1.d0
      go to 130
  200 lx=l-1
      sexit=1.d0
  210 do 310 i=1,3
      sum2=1.d0
      sum1=1.d0
      term1=1.d0
      term2=1.d0
      fm=1.d0
      tl=2*lx
  220 tm=2.d0*fm
      term1=-term1*zsq/(tm*(tl+tm+1.d0))*sexit
      term2=term2*zsq/(tm*(tl-tm+1.d0))*sexit
      old1=sum1
      old2=sum2
      sum1=sum1+term1
      sum2=sum2+term2
      diff=abs((old1-sum1)/sum1)+abs((old2-sum2)/sum2)
      if (diff-.000001d0) 240,240,230
  230 fm=fm+1.d0
      go to 220
  240 ix=2*lx-1
      asm=-1.d0
      if (lx) 270,270,250
  250 do 260 j=1,ix,2
      aj=j
  260 asm=asm*aj
  270 zyl=(z/y)**lx
      zyl1=zyl*(z/y)
      if (lx) 280,280,290
  280 asm1=1.d0
      go to 300
  290 alx=lx
      asm1=-asm*(2.d0*alx+1.d0)
  300 bj(i)=zyl*sum1/asm1
      by(i)=asm*sum2/zyl1
  310 lx=lx+1
      if (sexit) 340,340,570
  320 sign=-1.d0
      er=exp(zb)
      erl=1.d0/er
      snh=.5d0*(er-erl)
      csh=.5d0*(er+erl)
      bj(1)=snh/zb
      by(1)=y*csh
      go to 100
  330 lx=l-1
      sexit=-1.d0
      go to 210
  340 knt=2
      by(1)=y/zb
      by(2)=y*y*(1.d0+zb)/(zb*zb)
  350 by(3)=(2.d0*float(knt)-1.d0)*y*by(2)/zb+y*y*by(1)
      if (knt-l-1) 360,370,370
  360 by(1)=by(2)
      by(2)=by(3)
      knt=knt+1
      go to 350
  370 exz=exp(-zb)
      do 380 n=1,3
  380 by(n)=by(n)*exz
      if (zb-1.d0) 400,390,390
  390 if (z-4.d0*flr) 500,400,400
  400 er=2.d0*fl+1.d0
      bjp=(fl*bj(1)/y+fl1*bj(3)*y)/er
      byp=-(fl*by(1)*y+fl1*by(3)/y)/er
      bjpp=2.d0*bjp/z-(1.d0+fl*fl1/zsq)*bj(2)
      bypp=2.d0*byp/z-(1.d0+fl*fl1/zsq)*by(2)
      go to 580
  410 a=0.1d0
      ize=0
      b=0.35
      l1=l+1
      fn1=l
      idx=-1
      fn2=z-0.5d0+sqrt(30.0d0*b*z)
      ze=0.d0
      u1=2.d0*z/(2.d0*fn1+1.d0)
      u2=2.d0*z/(2.d0*fn2+1.d0)
      m1=fn1+30.d0*(a+b*u1*(2.d0-u1*u1)/(2.d0*(1.d0-u1*u1)))
      m2=fn2+30.d0*(a+b*u2*(2.d0-u2*u2)/(2.d0*(1.d0-u2*u2)))
      al=l
      if (fn2-al) 420,420,430
  420 am=m1+1
      m=m1
      go to 450
  430 if (m1-m2) 420,420,440
  440 am=m2+1
      m=m2
  450 r(3)=z/(2.d0*am+3.d0)
      n=m
  460 an=n
      r(2)=z/(2.d0*an+3.d0-z*r(3))
      n=n-1
      if (l-n) 470,470,480
  470 r(3)=r(2)
      go to 460
  480 bj(3)=r(2)
      bj(2)=1.d0
      al=2*l+1
      bj(1)=al/z-r(2)
      la=-l
      la1=-l-1
      alpha=z*z*(bj(2)*by(1)*y**la-bj(1)*by(2)*y**la1)
      na=l-1
      do 490 n=1,3
      bj(n)=(1.d0/(y**na*alpha))*bj(n)
  490 na=na+1
      go to 570
  500 a=0.1
      b=.35
      fn=-.5d0+sqrt(30.d0*b*zb)
      if (fn-fl) 510,520,520
  510 fn=fl
  520 u=2.d0*zb/(2.d0*fn+1.d0)
      m=fn+30.d0*(a+b*u)+1.d0
      r(3)=zb/2.d0*float(m)+1.d0
      m=m-1
  530 an=m
      r(2)=zb/(2.d0*an+1.d0+zb*r(3))
      m=m-1
      if (m-l-1) 550,540,540
  540 r(3)=r(2)
      go to 530
  550 bj(3)=r(2)
      bj(2)=1.d0
      bj(1)=(2.d0*fl+1.d0)/zb+r(2)
      alpha=zb*zb*(bj(2)*by(1)/y**l+bj(1)*by(2)/y**l1)
      na=fl-1.d0
      do 560 n=1,3
      bj(n)=(1.d0/(y**na*alpha))*bj(n)
  560 na=na+1
      go to 400
  570 er=1.d0/(2.d0*fl+1.d0)
      bjp=er*(fl*bj(1)/y-fl1*bj(3)*y)
      byp=er*(fl*by(1)*y-fl1*by(3)/y)
      bjpp=(fl*fl1/zsq-1.d0)*bj(2)-2.d0*bjp/z
      bypp=(fl*fl1/zsq-1.d0)*by(2)-2.d0*byp/z
  580 ba=bj(2)
      bb=by(2)
      y=ysave
      return
      end
