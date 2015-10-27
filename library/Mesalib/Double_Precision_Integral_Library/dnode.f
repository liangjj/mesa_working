*deck @(#)dnode.f	5.1 11/6/94 
      subroutine dnode(a,rt,k)
c
c
c     returns in rt(i) the ith root of a polynomial of order k whose
c     mth coefficient is stored in a(m+1). it is assumed that the
c     initial values in rt bracket the final values.
c
c
      implicit double precision (a-h,o-z)
cvax  implicit double precision (a-h,o-z)
c
      integer wptoin
      dimension a(10),rt(10)
c
      data a0,tm21,a1s16,a1s4,a3s4 /0.0d+00,1.0d-21,6.25d-02,0.25d+00,
     #                             0.75d+00/
      save a0,tm21,a1s16,a1s4,a3s4
c..bhl
      tol=tm21
c..bhl
      if(wptoin(1).eq.2) tol=1.0d-13
c..bhl
c.temp      tol=1.d-18
c..bhl
      k1=k+1
      r2=a0
      p2=a(1)
      do 100 m=1,k
      r1=r2
      p1=p2
      r2=rt(m)
      p2=a(k1)
      do 10 i=1,k
   10 p2=p2*r2+a(k1-i)
      prod=p1*p2
      if(prod.lt.a0) go to 20
      write(6,15) m,k
   15 format(/12h0root number,i4,
     1   38h was not found for polynomial of order,i4//)
      call lnkerr('stop in dnode')
   20 r5=r1
      p5=p1
      r6=r2
      p6=p2
   30 r3=r5
      p3=p5
      r4=r6
      p4=p6
      r =(r3*p4-r4*p3)/(p4-p3)
      dr=r4-r3
      delta=dr
      if(abs(delta).lt.tol) go to 90
      dr=a1s16*dr
      r5=r-dr
      if(r5.lt.r3) r5=r3
      r6=r+dr
      if(r6.gt.r4) r6=r4
      p5=a(k1)
      p6=p5
      do 40 i=1,k
      p5=p5*r5+a(k1-i)
   40 p6=p6*r6+a(k1-i)
   45 prod=p5*p6
      if(prod.lt.a0) go to 30
      prod=p3*p5
      if(prod.gt.a0) go to 60
      r5=a1s4*r3+a3s4*r5
      p5=a(k1)
      do 50 i=1,k
   50 p5=p5*r5+a(k1-i)
      go to 45
   60 r6=a1s4*r4+a3s4*r6
      p6=a(k1)
      do 70 i=1,k
   70 p6=p6*r6+a(k1-i)
      go to 45
   90 rt(m)=r
  100 continue
      return
      end
