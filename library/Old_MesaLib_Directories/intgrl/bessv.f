*deck @(#)bessv.f	5.1  11/6/94
      subroutine bessv(b,z,l,n)
      implicit real*8(a-h,o-z)
c
c     this routine evaluates modified spherical bessel functions.
c
c     the arguments:
c       b...  bessel function at requested points.
c       z...  the points at which the function is to be evaluated.
c             assumed to be in ascending order.
c       l...  the l value of the function.
c       n...  the number of requested points.
c
c
c       for arguments (z.le.5.0)  use power series
c                                 15 terms good to 5.0d-14
c                     (z.gt.5.0)  use exponential representation.
c                     (z.gt.16.1) only first term in exponential
c                                 representation is required.
c
c     result has a factor of exp(-z) included to avoid overflow.
c
      common/dfac/dfac(30)
      common/fact/fact(17),fprod(9,9)
      common/const/zero,one,two,three,four,five,six,ten
c
      real*8 z(1),b(1)
      real*8 den(15),term(50),zp(50)
      real*8 rp(50),rm(50),tzpi(50),tzmi(50),texm(50)
      real*8 denm(50,9),denp(50,9)
c
      equivalence (den(1),rp(1)), (term(1),rm(1)), (zp(1),tzpi(1))
c
      data pt5/0.5d0/,f16pt1/16.1d0/
      save pt5,f16pt1
c
c     ----- determine the number of points in each region -----
      ilo1=0
      ilo2=0
      ilo3=0
      n1=0
      n2=0
      n3=0
      do 10 i=1,n
         if(z(i).lt.zero) then
            b(i)=zero
         else if (z(i).eq.zero) then
            b(i)=zero
            if(l.eq.0) b(i)=one
         else
            if(z(i).le.five) then
c              0<= z <5.0
               if(ilo1.eq.0) ilo1=i
               n1=n1+1
            else if (z(i).le.f16pt1) then
c              5.0< z <16.1
               if(ilo2.eq.0) ilo2=i
               n2=n2+1
            else
c              16.1 < z
               if(ilo3.eq.0) ilo3=i
               n3=n3+1
            endif
         endif
   10 continue
      ihi1=ilo1+n1-1
      ihi2=ilo2+n2-1
      ihi3=ilo3+n3-1
c
      if(n1.ne.0) then
c        power series.
         do 20 i=ilo1,ihi1
            zp(i)=z(i)*z(i)*pt5
            term(i)=(z(i)**l)/dfac(l+l+3)
            b(i)=term(i)
   20    continue
         do 23 j=1,15
            den(j)=one/float(j*(l+l+j+j+1))
   23    continue
         do 25 j=1,15
            do 25 i=ilo1,ihi1
               term(i)=term(i)*zp(i)*den(j)
               b(i)=b(i)+term(i)
   25    continue
         do 28 i=ilo1,ihi1
            b(i)=b(i)*exp(-z(i))
   28    continue
      endif
c
c     exponential represention.
      if(n2.ne.0) then
         l1=l+1
         do 52 i=ilo2,ihi2
            rp(i)=fprod(1,l1)
            rm(i)=fprod(1,l1)
            tzi=two*z(i)
            tzpi(i)=one/tzi
            tzmi(i)=-tzpi(i)
            texm(i)=exp(-tzi)
            denp(i,1)=one
            denm(i,1)=one
   52    continue
         if(l1.gt.1) then
            do 54 k1=2,l1
               do 54 i=ilo2,ihi2
                  denp(i,k1)=denp(i,k1-1)*tzpi(i)
                  denm(i,k1)=denm(i,k1-1)*tzmi(i)
   54       continue
            do 60 k1=2,l1
               do 60 i=ilo2,ihi2
                  rp(i)=rp(i)+denp(i,k1)*fprod(k1,l1)
                  rm(i)=rm(i)+denm(i,k1)*fprod(k1,l1)
   60       continue
         endif
         do 62 i=ilo2,ihi2
            b(i)=(rm(i)-((-one)**l)*rp(i)*texm(i))*tzpi(i)
   62    continue
      endif
c
c     only first term in exponential representation.
      if(n3.ne.0) then
         l1=l+1
         do 72 i=ilo3,ihi3
            rm(i)=fprod(1,l1)
            tzmi(i)=-pt5/z(i)
            denm(i,1)=one
   72    continue
         if(l1.gt.1) then
            do 74 k1=2,l1
               do 74 i=ilo3,ihi3
                  denm(i,k1)=denm(i,k1-1)*tzmi(i)
   74       continue
            do 80 k1=2,l1
               do 80 i=ilo3,ihi3
                  rm(i)=rm(i)+denm(i,k1)*fprod(k1,l1)
   80       continue
         endif
         do 82 i=ilo3,ihi3
            b(i)=-rm(i)*tzmi(i)
   82    continue
      endif
c
c
      return
      end
