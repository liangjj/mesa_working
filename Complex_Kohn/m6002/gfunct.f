*deck @(#)gfunct.f	1.1 9/7/91
c***begin prologue     gfunct
c***date written                (yymmdd)
c***revision date      890411   (yymmdd)
c***keywords           gfunct, link 6003
c***authors            unknown
c***                   
c***source             m6002
c***purpose            special functions for static potential
c***references       
c
c***routines called    
c***end prologue       gfunct
      subroutine gfunct (l,m,a,b,p,t,g,n,narg)
      implicit real *8 (a-h,o-z)
      dimension p(narg), g(narg,7,3)
      ll=l+1
      mm=m+1
      if (ll.eq.1) then
          if (mm.eq.1) then
              do 210 ig=1,narg
                 g(ig,1,n)=1.d0
  210         continue
          elseif (mm.eq.2) then
              do 211 ig=1,narg
                 g(ig,1,n)=b
                 g(ig,2,n)=-p(ig)
  211         continue
          elseif (mm.eq.3) then
              do 212 ig=1,narg
                 g(ig,1,n)=b*b+0.5d0*t
                 g(ig,2,n)=-2.d0*b*p(ig)-0.5d0*t
                 g(ig,3,n)=p(ig)*p(ig)
  212         continue 
          elseif (mm.eq.4) then
              temp=a
              a=b
              b=temp
              do 240 i=1,narg
                 g(i,1,n)=a*(a*a+1.5d0*t)
                 g(i,2,n)=-3.d0*(a*(a*p(i)+0.5d0*t)+0.5d0*p(i)*t)
                 g(i,3,n)=3.d0*p(i)*(a*p(i)+0.5d0*t)
                 g(i,4,n)=-p(i)*p(i)*p(i)
240           continue
           else
              call lnkerr('error in gfunct routine. ll=1')
           endif 
      elseif (ll.eq.2) then
          if (mm.eq.1) then         
              do 220 ig=1,narg
                 g(ig,1,n)=a
                 g(ig,2,n)=-p(ig)
220           continue
          elseif (mm.eq.2) then
              do 221 ig=1,narg
                 g(ig,1,n)=a*b+0.5d0*t
                 g(ig,2,n)=-p(ig)*(a+b)-0.5d0*t
                 g(ig,3,n)=p(ig)*p(ig)
221           continue
          elseif (mm.eq.3) then
              do 222 ig=1,narg
                 g(ig,1,n)=b*b*a+0.5d0*t*(a+2.d0*b)
                 g(ig,2,n)=-p(ig)*b*(2.d0*a+b)-0.5d0*t*((a+2.d0*b)+
     1                      3.d0*p(ig))
                 g(ig,3,n)=p(ig)*(p(ig)*(a+2.d0*b)+1.5d0*t)
                 g(ig,4,n)=-p(ig)*p(ig)*p(ig)
222           continue
          elseif(mm.eq.4) then
              temp=a
              a=b
              b=temp
              t2=t*t
              a2=a*a
              ab=a*b
              f0=a2*ab
              f1=a2*(a+3.d0*b)
              f2=3.d0*a*(a+b)
              f3=3.d0*a+b
              do 241 i=1,narg
                 p2=p(i)*p(i)
                 g(i,1,n)=f0+.5d0*t*f2+.75d0*t2
                 g(i,2,n)=-(p(i)*f1+.5d0*t*f2+1.5d0*p(i)*t*f3+1.5d0*t2)
                 g(i,3,n)=p2*f2+1.5d0*p(i)*t*f3+3.d0*t*p2+.75d0*t2
                 g(i,4,n)=-(p(i)*p2*f3+3.d0*t*p2)
241           continue
          else
             call lnkerr('error in gfunct routine. ll=2')
          endif 
      elseif(ll.eq.3) then
          if (mm.eq.1) then
              do 230 ig=1,narg
                 g(ig,1,n)=a*a+0.5d0*t
                 g(ig,2,n)=-2.d0*a*p(ig)-0.5d0*t
                 g(ig,3,n)=p(ig)*p(ig)
230           continue
          elseif (mm.eq.2) then
              do 231 ig=1,narg
                 g(ig,1,n)=a*a*b+0.5d0*t*(2.d0*a+b)
                 g(ig,2,n)=-p(ig)*a*(a+2.d0*b)-0.5d0*t*((2.d0*a+b)+
     1                      3.d0*p(ig))
                 g(ig,3,n)=p(ig)*(p(ig)*(2.d0*a+b)+1.5d0*t)
                 g(ig,4,n)=-p(ig)*p(ig)*p(ig)
231           continue
          elseif (mm.eq.3) then
              aa=a*a
              bb=b*b
              ab=4.d0*a*b
              do 232 i=1,narg
                 pp=p(i)*p(i)
                 g(i,1,n)=aa*bb+t*(0.5d0*(aa+ab+bb)+0.75d0*t)
                 g(i,2,n)=-(2.d0*p(i)*(aa*b+a*bb)+t*(0.5d0*(aa+ab+bb)+
     1                      3.d0*(a+b)*p(i)+1.5d0*t))
                 g(i,3,n)=pp*((aa+ab+bb)+3.d0*t)+t*(3.d0*p(i)*(a+b)+
     1                        0.75d0*t)
                 g(i,4,n)=-(pp*(2.d0*p(i)*(a+b)+3.d0*t))
                 g(i,5,n)=pp*pp
232           continue
          elseif (mm.eq.4) then
              temp=a
              a=b
              b=temp
              a2=a*a
              b2=b*b
              ab=a*b
              t2=t*t
              a3=3.d0*a2+6.d0*ab+b2
              a1=a2+6.d0*ab+3.d0*b2
              do 242 i=1,narg
                 p2=p(i)*p(i)
                 g(i,1,n)=a*a2*b2+0.5d0*t*a1*a+0.75d0*t2*(3.d0*a+
     1                    2.d0*b)
                 g(i,2,n)=-(ab*p(i)*(2.d0*a2+3.d0*ab)+0.5d0*t*a*a1+
     1                      1.5d0*t*a3*p(i)+1.5d0*t2*(3.d0*a+2.d0*b)+
     2                      3.75d0*p(i)*t2)
                 g(i,3,n)=p2*a*a1+1.5d0*p(i)*t*a3+3.d0*(3.d0*a+2.d0*b)*
     1                    (t*p2+.25d0*t2)+7.5d0*p(i)*t2
                 g(i,4,n)=-(p2*(p(i)*a3+t*(9.d0*a+6.d0*b)+5.d0*p(i)*t)+
     1                      3.75d0*p(i)*t2)
                 g(i,5,n)=p(i)*p2*(p(i)*(3.d0*a+2.d0*b)+5.d0*t)
                 g(i,6,n)=-p(i)*p2*p2
242           continue
          else
             call lnkerr('error in gfunct routine. ll=3')
          endif 
      elseif(ll.eq.4) then
          if (mm.eq.1) then
              do 340 i=1,narg
                 g(i,1,n)=a*(a*a+1.5d0*t)
                 g(i,2,n)=-3.d0*(a*(a*p(i)+0.5d0*t)+0.5d0*p(i)*t)
                 g(i,3,n)=3.d0*p(i)*(a*p(i)+0.5d0*t)
                 g(i,4,n)=-p(i)*p(i)*p(i)
340           continue
          elseif (mm.eq.2) then
              t2=t*t
              a2=a*a
              ab=a*b
              f0=a2*ab
              f1=a2*(a+3.d0*b)
              f2=3.d0*a*(a+b)
              f3=3.d0*a+b
              do 341 i=1,narg
                 p2=p(i)*p(i)
                 g(i,1,n)=f0+.5d0*t*f2+.75d0*t2
                 g(i,2,n)=-(p(i)*f1+.5d0*t*f2+1.5d0*p(i)*t*f3+1.5d0*t2)
                 g(i,3,n)=p2*f2+1.5d0*p(i)*t*f3+3.d0*t*p2+.75d0*t2
                 g(i,4,n)=-(p(i)*p2*f3+3.d0*t*p2)
                 g(i,5,n)=p2*p2      
341           continue
          elseif (mm.eq.3) then 
              a2=a*a
              b2=b*b
              ab=a*b
              t2=t*t
              a3=3.d0*a2+6.d0*ab+b2
              a1=a2+6.d0*ab+3.d0*b2
              do 342 i=1,narg
                 p2=p(i)*p(i)
                 g(i,1,n)=a*a2*b2+0.5d0*t*a1*a+0.75d0*t2*(3.d0*a+
     1                    2.d0*b)
                 g(i,2,n)=-(ab*p(i)*(2.d0*a2+3.d0*ab)+0.5d0*t*a*a1+
     1                      1.5d0*t*a3*p(i)+1.5d0*t2*(3.d0*a+2.d0*b)+
     2                      3.75d0*p(i)*t2)
                 g(i,3,n)=p2*a*a1+1.5d0*p(i)*t*a3+3.d0*(3.d0*a+2.d0*b)*
     1                    (t*p2+.25d0*t2)+7.5d0*p(i)*t2
                 g(i,4,n)=-(p2*(p(i)*a3+t*(9.d0*a+6.d0*b)+5.d0*p(i)*t)+
     1                      3.75d0*p(i)*t2)
                 g(i,5,n)=p(i)*p2*(p(i)*(3.d0*a+2.d0*b)+5.d0*t)
                 g(i,6,n)=-p(i)*p2*p2
342           continue
          elseif (mm.eq.4) then
              a2=a*a
              b2=b*b
              t2=t*t
              ab=a*b
              f0=a2*b2*ab
              f1=3.d0*a2*b2*(a+b)
              f2=3.d0*ab*(a2+3.d0*ab+b2)
              f3=a2*(a+9.d0*b)+b2*(b+9.d0*a)
              f4=3.d0*(a2+3.d0*ab+b2)
              f5=3.d0*(a+b)
              do 343 i=1,narg
                 p2=p(i)*p(i)
                 g(i,1,n)=f0+.5d0*t*f2+.75d0*t2*f4+1.875d0*t*t2
                 g(i,2,n)=-(p(i)*f1+.5d0*t*f2+1.5d0*t*p(i)*f3+
     1                      1.5d0*t2*f4+3.75d0*p(i)*t2*f5+
     2                      5.625d0*t*t2)
                 g(i,3,n)=p2*f2+1.5d0*p(i)*t*f3+3.d0*f4*(t*p2+
     1                   .25d0*t2)+7.5d0*p(i)*t2*f5+11.25d0*t2*
     2                    (.5d0*t+p2)
                 g(i,4,n)=-(p2*(p(i)*f3+3.d0*t*f4)+5.d0*f5*p(i)*t*
     1                     (p2+.75d0*t)+t2*(22.5d0*p2+1.875d0*t))
                 g(i,5,n)=p2*(f4*p2+5.d0*p(i)*t*f5+7.5d0*p2*t+
     1                    11.25d0*t2)
                 g(i,6,n)=-(p2*p2*(p(i)*f5+7.5d0*t))
                 g(i,7,n)=p2*p2*p2
343           continue
          else
             call lnkerr('error in gfunct routine. ll=4')
          endif 
      else
          call lnkerr('error in gfunct routine. ll out of range')
      endif
      return
      end 
