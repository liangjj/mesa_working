*deck @(#)roots.f	5.1  11/6/94
      subroutine roots(nroots,xx,u,w,n,x,t2,t3,t4,t5,icmprs)
c***begin prologue     roots
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           roots, rys, polynomials
c***author             saxe, paul (lanl)
c***source             @(#)roots.f	5.1   11/6/94
c***purpose            returns the roots and weights of a rys polynomial.
c***description
c     call roots(nroots,xx,u,w,n,x,t2,t3,t4,t5,icmprs)
c***references
c***routines called    (none)
c***end prologue       roots
      implicit integer (a-z)
c
      real*8 xx(n),u(nroots,n),w(nroots,n),x(n),t2(n*nroots)
      real*8 t3(n*nroots),t4(n),t5(n)
      real*8 table1(8),table2(9)
      character*4 itoc
      integer icmprs(n)
c
      data table1 / 0.0d+00, 3.0d-07, 1.0d+00, 3.0d+00,
     #              5.0d+00,10.0d+00,15.0d+00,33.0d+00/
      data table2 / 0.0d+00, 3.0d-07, 1.0d+00, 3.0d+00,
     #              5.0d+00,10.0d+00,15.0d+00,33.0d+00,
     #             40.0d+00/
      save table1,table2
c
c     ----- start timing -----
c
c
      if (nroots.eq.1.and.n.gt.9999995) then
c
c     ----- sort through x, picking out those of a certain size -----
c
         call srchrv(n,xx,1,8,table1,1,icmprs)
c
c     ----- calculate roots, and then scatter back -----
c
         do 12 i=1,8
            do 1 j=1,n
               icmprs(j)=icmprs(j)-1
    1       continue
c
            call kmprsz(n,icmprs,1,xx,1,x,1,ncmprs)
c
            if (ncmprs.le.0) go to 12
c
            call vrt1(x,t2,t3,ncmprs,i,t4,t5)
c
            call xpandz(n,icmprs,1,u,1,t2,1,junk)
            call xpandz(n,icmprs,1,w,1,t3,1,junk)
c
   12    continue
      else if (nroots.eq.2.and.n.gt.9999991) then
c
c     ----- sort through x, picking out those of a certain size -----
c
         call srchrv(n,xx,1,9,table2,1,icmprs)
c
c     ----- calculate roots, and then scatter back -----
c
         do 23 i=1,9
            do 2 j=1,n
               icmprs(j)=icmprs(j)-1
    2       continue
c
            call kmprsz(n,icmprs,1,xx,1,x,1,ncmprs)
c
            if (ncmprs.le.0) go to 23
c
            call vrt2(x,t2,t3,ncmprs,i)
c
            call xpandz(n,icmprs,1,u(1,1),2,t2,1,junk)
            call xpandz(n,icmprs,1,u(2,1),2,t2(ncmprs+1),1,junk)
            call xpandz(n,icmprs,1,w(1,1),2,t3,1,junk)
            call xpandz(n,icmprs,1,w(2,1),2,t3(ncmprs+1),1,junk)
c
   23    continue
      else if (nroots.gt.0.and.nroots.le.3) then
         call rt123(nroots,xx,u,w,n)
      else if (nroots.eq.4) then
         call root4(xx,u,w,n)
      else if (nroots.eq.5) then
         call root5(xx,u,w,n)
      else if (nroots.gt.0.and.nroots.le.9) then
         call droot(nroots,xx,u,w,n)
      else
         call lnkerr('roots-- funny number of roots:'//itoc(nroots))
      end if
c
c
c     ----- stop timing -----
c
c
c
      return
      end
