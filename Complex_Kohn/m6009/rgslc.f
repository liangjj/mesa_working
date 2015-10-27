*deck @(#)rgslc.f	1.1 9/8/91
c***begin prologue     rgslc
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           rgslc, link m1106, kohn variational
c***author             schneider, barry (lanl), rescigno, tom(llnl)  
c***source             m1106
c***purpose            solve complex linear systems
c***description        solve a linear equation in which
c***                   the matrix has been factored, is real
c***                   and the right hand side is complex.
c***                   this is a specialization of the linpack
c***                   routine sgesl to a complex rhs.                  
c***references         schneider and rescigno, physical review
c
c***routines called    iosys, util and mdutil
c***end prologue       rgslc
      subroutine rgslc(a,lda,n,ipvt,b)
      implicit integer (a-z)
      real *8 a
      complex*16  b, t
      dimension a(lda,n), b(n), ipvt(n)
      nm1=n-1
      if (nm1.lt.1) go to 30
          do 20 k=1,nm1
             l=ipvt(k)
             t=b(l)
             if (l.eq.k) go to 10
                 b(l)=b(k)
                 b(k)=t
   10        continue
             call csxpy(n-k,t,a(k+1,k),b(k+1))
   20     continue
   30 continue
         do 40 kb=1,n
             k=n+1-kb
             b(k)=b(k)/a(k,k)
             t=-b(k)
             call csxpy(k-1,t,a(1,k),b(1))
   40     continue
      return
      end
