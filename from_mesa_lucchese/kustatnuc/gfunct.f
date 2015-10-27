      subroutine gfunct (l,m,a,b,p,t,g,n,narg)
      parameter (mxbuf=1000)
      implicit real*8 (a-h,o-z)
      dimension p(mxbuf),g(mxbuf,7,3)
      ll=l+1
      mm=m+1
      go to (100,101,102,103),ll
  100 go to (110,111,112,113),mm
  101 go to (120,121,122,123),mm
  102 go to (130,131,132,133),mm
  103 go to (140,141,142,143),mm
110   do 210 ig=1,narg
210   g(ig,1,n)=1.e0
      go to 300
111   do 211 ig=1,narg
      g(ig,1,n)=b
  211 g(ig,2,n)=-p(ig)
      go to 300
  112 do 212 ig=1,narg
      g(ig,1,n)=b*b+0.5e0*t
      g(ig,2,n)=-2.e0*b*p(ig)-0.5e0*t
  212 g(ig,3,n)=p(ig)*p(ig)
      go to 300
  113 temp=a
      a=b
      b=temp
      go to 140
120   do 220 ig=1,narg
      g(ig,1,n)=a
      g(ig,2,n)=-p(ig)
220   continue
      go to 300
121   do 221 ig=1,narg
      g(ig,1,n)=a*b+0.5e0*t
      g(ig,2,n)=-p(ig)*(a+b)-0.5e0*t
      g(ig,3,n)=p(ig)*p(ig)
221   continue
      go to 300
122   do 222 ig=1,narg
      g(ig,1,n)=b*b*a+0.5e0*t*(a+2.e0*b)
      g(ig,2,n)=-p(ig)*b*(2.e0*a+b)-0.5e0*t*((a+2.e0*b)+3.e0*p(ig))
      g(ig,3,n)=p(ig)*(p(ig)*(a+2.e0*b)+1.5e0*t)
      g(ig,4,n)=-p(ig)*p(ig)*p(ig)
222   continue
      go to 300
  123 temp=a
      a=b
      b=temp
      go to 141
130   do 230 ig=1,narg
      g(ig,1,n)=a*a+0.5e0*t
      g(ig,2,n)=-2.e0*a*p(ig)-0.5e0*t
      g(ig,3,n)=p(ig)*p(ig)
230   continue
      go to 300
131   do 231 ig=1,narg
      g(ig,1,n)=a*a*b+0.5e0*t*(2.e0*a+b)
      g(ig,2,n)=-p(ig)*a*(a+2.e0*b)-0.5e0*t*((2.e0*a+b)+3.e0*p(ig))
      g(ig,3,n)=p(ig)*(p(ig)*(2.e0*a+b)+1.5e0*t)
      g(ig,4,n)=-p(ig)*p(ig)*p(ig)
231   continue
      go to 300
  132 aa=a*a
      bb=b*b
      ab=4.e0*a*b
      do 232 i=1,narg
      pp=p(i)*p(i)
      g(i,1,n)=aa*bb+t*(0.5e0*(aa+ab+bb)+0.75e0*t)
      g(i,2,n)=-(2.e0*p(i)*(aa*b+a*bb)+t*(0.5e0*(aa+ab+bb)+3.e0*(a+b)
     x*p(i)+1.5e0*t))
      g(i,3,n)=pp*((aa+ab+bb)+3.e0*t)+t*(3.e0*p(i)*(a+b)+0.75e0*t)
      g(i,4,n)=-(pp*(2.e0*p(i)*(a+b)+3.e0*t))
      g(i,5,n)=pp*pp
232   continue
      go to 300
  133 temp=a
      a=b
      b=temp
      go to 142
140   do 240 i=1,narg
      g(i,1,n)=a*(a*a+1.5e0*t)
      g(i,2,n)=-3.e0*(a*(a*p(i)+0.5e0*t)+0.5e0*p(i)*t)
      g(i,3,n)=3.e0*p(i)*(a*p(i)+0.5e0*t)
      g(i,4,n)=-p(i)*p(i)*p(i)
240   continue
      go to 300
141   continue
      t2=t*t
      a2=a*a
      ab=a*b
      f0=a2*ab
      f1=a2*(a+3.e0*b)
      f2=3.e0*a*(a+b)
      f3=3.e0*a+b
      do 241 i=1,narg
      p2=p(i)*p(i)
      g(i,1,n)=f0+.5e0*t*f2+.75e0*t2
      g(i,2,n)=-(p(i)*f1+.5e0*t*f2+1.5e0*p(i)*t*f3+1.5e0*t2)
      g(i,3,n)=p2*f2+1.5e0*p(i)*t*f3+3.e0*t*p2+.75e0*t2
      g(i,4,n)=-(p(i)*p2*f3+3.e0*t*p2)
c..ohio fix for f-functions
      g(i,5,n)=p2*p2
c..ohio fix for f-functions
241   continue
      go to 300
142   continue
      a2=a*a
      b2=b*b
      ab=a*b
      t2=t*t
      a3=3.e0*a2+6.e0*ab+b2
      a1=a2+6.e0*ab+3.e0*b2
      do 242 i=1,narg
      p2=p(i)*p(i)
      g(i,1,n)=a*a2*b2+0.5e0*t*a1*a+0.75e0*t2*(3.e0*a+2.e0*b)
      g(i,2,n)=-(ab*p(i)*(2.e0*a2+3.e0*ab)+0.5e0*t*a*a1+1.5e0*t*a3*p(i)+
     x 1.5e0*t2
     x*(3.e0*a+2.e0*b)+3.75e0*p(i)*t2)
      g(i,3,n)=p2*a*a1+1.5e0*p(i)*t*a3+3.e0*(3.e0*a+2.e0*b)*(t*p2+.25e0*
     x t2)+7
     x.5e0*p(i)*t2
      g(i,4,n)=-(p2*(p(i)*a3+t*(9.e0*a+6.e0*b)+5.e0*p(i)*t)+3.75e0*p(i)*
     x t2)
      g(i,5,n)=p(i)*p2*(p(i)*(3.e0*a+2.e0*b)+5.e0*t)
      g(i,6,n)=-p(i)*p2*p2
242   continue
      go to 300
143   continue
      a2=a*a
      b2=b*b
      t2=t*t
      ab=a*b
      f0=a2*b2*ab
      f1=3.e0*a2*b2*(a+b)
      f2=3.e0*ab*(a2+3.e0*ab+b2)
      f3=a2*(a+9.e0*b)+b2*(b+9.e0*a)
      f4=3.e0*(a2+3.e0*ab+b2)
      f5=3.e0*(a+b)
      do 243 i=1,narg
      p2=p(i)*p(i)
      g(i,1,n)=f0+.5e0*t*f2+.75e0*t2*f4+1.875e0*t*t2
      g(i,2,n)=-(p(i)*f1+.5e0*t*f2+1.5e0*t*p(i)*f3+1.5e0*t2*f4+3.75e0
     x*p(i)*t2*f5+5.
     x625e0*t*  t2)
      g(i,3,n)=p2*f2+1.5e0*p(i)*t*f3+3.e0*f4*(t*p2+.25e0*t2)+7.5e0*p(i)*
     x t2*f5+11
     x.25e0*t2* (.5e0*t+p2)
      g(i,4,n)=-(p2*(p(i)*f3+3.e0*t*f4)+5.e0*f5*p(i)*t*(p2+.75e0*t)+t2*(
     x 22.5e0*
     xp2+1.875e0*t))
      g(i,5,n)=p2*(f4*p2+5.e0*p(i)*t*f5+7.5e0*p2*t+11.25e0*t2)
      g(i,6,n)=-(p2*p2*(p(i)*f5+7.5e0*t))
      g(i,7,n)=p2*p2*p2
243   continue
  300 continue
      return
      end
