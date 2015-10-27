      subroutine pts(nalpha,nbeta,ngamma,x,w)
      implicit real*8(a-h,o-z)
      real*8 nmula,nmulb,nmulg
      dimension xx(64),ww(64),x(64,3),w(64,3)
     1,dummy(2),scr(64)
      read(5,*)itype,jtype,ktype
      read(5,*)nalpha,nbeta,ngamma
      read(5,*)a1,a2
      read(5,*)b1,b2
      read(5,*)g1,g2
      pi=3.14159265358
      fac=pi/180.
      a1=fac*a1
      a2=fac*a2
      b1=fac*b1
      b2=fac*b2
      g1=fac*g1
      g2=fac*g2
      faca=a2-a1
      facb=b2-b1
      facg=g2-g1
      nmula=2.*pi/faca
      nmulb=pi/facb
      nmulg=2*pi/facg
      write(6,16)nmula,nmulb,nmulg
16    format(" multiplicative factors for euler angles  :",3f10.5)
      irun=0
      call gaussq(itype,nalpha,0.,0.,0,dummy,scr,xx,ww)
      ampa=(a2-a1)/2.
      abpa=(a2+a1)/2.
      do 1 i=1,nalpha
      x(i,1)=ampa*xx(i)+abpa
1     w(i,1)=nmula*ampa*ww(i)
      call gaussq(jtype,nbeta,0.,0.,0,dummy,scr,xx,ww)
      ampb=(cos(b1)-cos(b2))/2.
      abpb=(cos(b2)+cos(b1))/2.
      do 2 i=1,nbeta
      x(i,2)=ampb*xx(i)+abpb
2     w(i,2)=nmulb*ampb*ww(i)
      call gaussq(ktype,ngamma,0.,0.,0,dummy,scr,xx,ww)
      ampg=(g2-g1)/2.
      abpg=(g2+g1)/2.
      do 3 i=1,ngamma
      x(i,3)=ampg*xx(i)+abpg
3     w(i,3)=nmulg*ampg*ww(i)
      return
      end
