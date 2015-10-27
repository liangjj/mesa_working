*deck @(#)dfunc.f	5.1 11/6/94 
      subroutine dfunc(x,n)
c
c
      implicit double precision (a-h,o-z)
cvax  implicit double precision (a-h,o-z)
c
      integer wptoin
      common /ffm/    ff(19)
c
      data tm78, tm29, a1s2, pie4, a1, a2, a80, a360 /1.0d-78, 1.0d-29,
     1   0.5d+00,0.7853981633974483096156608d+00,1.0d+00,2.0d+00,
     #  80.0d+00,360.0d+00/
      data a0 /0.0d+00/
      save tm78,tm29,a1s2,pie4,a1,a2,a80,a360
      save a0
c
      tol=tm29
c
c     for 32-bit machines.
c
      if (wptoin(1).eq.2) tol=1.0d-14
c
      xx=x+x
      facmin=xx
      e=tm78
      qx=-x
      if(facmin.lt.a360) e=exp(qx)
      if (e.eq.a0) e=tm78
      if(facmin.gt.a80) go to 100
      term=a1
      sum=a1
      fac=n
      fac=fac+a1s2
   10 fac=fac+a1
      term=term*x/fac
      sum=sum+term
      if(fac.le.facmin) go to 10
      t=term
      s=sum
      if(t.gt.s*tol) go to 10
      fac=n+n+1
      ff(n+1)=sum*e/fac
      m=n-1
      fac=m+m+1
   20 if(m.lt.0) return
      ff(m+1)=(e+xx*ff(m+2))/fac
      m=m-1
      fac=fac-a2
      go to 20
c
c     use asymptotic expansion for large arguments.
c
  100 a=sqrt(pie4/x)
      tmax=a*tol/e
      term=a1/xx
      sum=term
      fac=a1
  110 fac=fac-a2
      term=fac*term/xx
      sum=term+sum
      t=term
      if(abs(t).gt.tmax) go to 110
      ff(1)=a-e*sum
      fac=-a1
      m=0
  120 if(m.eq.n) return
      m=m+1
      fac=fac+a2
      ff(m+1)=(fac*ff(m)-e)/xx
      go to 120
      end
