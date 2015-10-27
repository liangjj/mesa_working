      subroutine plm(x,n,mu,lmax,p,p0,p1,p2,arg,q0,q1,q2,ngrid)
      implicit real*8(a-h,o-z)
      dimension p(ngrid,0:lmax),x(ngrid),p0(ngrid),p1(ngrid),p2(ngrid)
      dimension arg(ngrid),q0(ngrid),q1(ngrid),q2(ngrid)
c
c this routine calculates a vector of associated legengre polynomials P(l,mu) .
c for fixed mu and l=mu,,,lmax for an array of n points, x(i),i=1,,n.
c
c
c mu = 0 case
c
      if(mu.eq.0)then
       do 2 i=1,n
       p(i,0)=1.
2      p(i,1)=x(i)
       if(lmax.lt.2)return
       do 1 l=2,lmax
       do 3 i=1,n
3      p(i,l)=((2*l-1)*x(i)*p(i,l-1)-(l-1)*p(i,l-2))/l
1      continue
       return
      endif
c
c mu = 1 case
c
      if(mu.eq.1)then
       do 6 i=1,n
       arg(i)=sqrt(1.-x(i)*x(i))
6      p(i,1)=-arg(i)
       if(lmax.lt.2)return
       do 7 i=1,n
7      p(i,2)=-3.*x(i)*arg(i)
       if(lmax.lt.3)return
       do 8 l=3,lmax
       do 8 i=1,n
8      p(i,l)=((2*l-1)*x(i)*p(i,l-1)-l*p(i,l-2))/(l-1)
       return
      endif
c
c mu must be larger than 1
c
c
c recurr across to l=mu+1, with m=0
c
      mu1=mu+1
      do 4 i=1,n
      p0(i)=1.
4     p1(i)=x(i)
      do 5 l=2,mu1
      do 5 i=1,n
      p2(i)=((2*l-1)*x(i)*p1(i)-(l-1)*p0(i))/l
      p0(i)=p1(i)
      p1(i)=p2(i)
5     continue
c
c recurr across to l=mu+1, with m=1
c
      do 9 i=1,n
      arg(i)=sqrt(1.-x(i)*x(i))
      q0(i)=-arg(i)
9     q1(i)=-3.*arg(i)*x(i)
      do 10 l=3,mu1
      do 10 i=1,n
      q2(i)=((2*l-1)*x(i)*q1(i)-l*q0(i))/(l-1)
      q0(i)=q1(i)
10    q1(i)=q2(i)
c
c with l fixed at mu and mu+1, recurr down to m=mu
c
      l0=mu
      l1=mu1
      do 11 m=2,mu
      do 11 i=1,n
      p2(i)=-2.*(m-1)*x(i)*q0(i)/arg(i)-(l0-m+2)*(l0+m-1)*p0(i)
      q2(i)=-2.*(m-1)*x(i)*q1(i)/arg(i)-(l1-m+2)*(l1+m-1)*p1(i)
      p0(i)=q0(i)
      p1(i)=q1(i)
      q0(i)=p2(i)
11    q1(i)=q2(i)
c
c compute the desired vector and quit
c
      do 12 i=1,n
      p(i,mu)=q0(i)
12    continue
      if(lmax.eq.mu)return
      do 13 i=1,n
13    p(i,mu1)=q1(i)
      if(lmax.eq.mu1)return
      mu2=mu1+1
      do 14 l=mu2,lmax
      do 14 i=1,n
14    p(i,l)=((2*l-1)*x(i)*p(i,l-1)-(l+mu-1)*p(i,l-2))/(l-mu)
      return
      end
