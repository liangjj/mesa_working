*deck hstcs1
      subroutine hstcs1 (intl, a, b, m, mbdcnd, bda, bdb, c, d, n,
     +   nbdcnd, bdc, bdd, elmbda, f, idimf, pertrb, ierr1, am, bm, cm,
     +   an, bn, cn, snth, rsq, wrk)
c***begin prologue  hstcs1
c***subsidiary
c***purpose  subsidiary to hstcsp
c***library   slatec
c***type      single precision (hstcs1-s)
c***author  (unknown)
c***see also  hstcsp
c***routines called  blktri
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  hstcs1
      dimension       bda(*)     ,bdb(*)     ,bdc(*)     ,bdd(*)     ,
     1                f(idimf,*) ,am(*)      ,bm(*)      ,cm(*)      ,
     2                an(*)      ,bn(*)      ,cn(*)      ,snth(*)    ,
     3                rsq(*)     ,wrk(*)
c***first executable statement  hstcs1
      dth = (b-a)/m
      dthsq = dth*dth
      do 101 i=1,m
         snth(i) = sin(a+(i-0.5)*dth)
  101 continue
      dr = (d-c)/n
      do 102 j=1,n
         rsq(j) = (c+(j-0.5)*dr)**2
  102 continue
c
c     multiply right side by r(j)**2
c
      do 104 j=1,n
         x = rsq(j)
         do 103 i=1,m
            f(i,j) = x*f(i,j)
  103    continue
  104 continue
c
c      define coefficients am,bm,cm
c
      x = 1./(2.*cos(dth/2.))
      do 105 i=2,m
         am(i) = (snth(i-1)+snth(i))*x
         cm(i-1) = am(i)
  105 continue
      am(1) = sin(a)
      cm(m) = sin(b)
      do 106 i=1,m
         x = 1./snth(i)
         y = x/dthsq
         am(i) = am(i)*y
         cm(i) = cm(i)*y
         bm(i) = elmbda*x*x-am(i)-cm(i)
  106 continue
c
c     define coefficients an,bn,cn
c
      x = c/dr
      do 107 j=1,n
         an(j) = (x+j-1)**2
         cn(j) = (x+j)**2
         bn(j) = -(an(j)+cn(j))
  107 continue
      isw = 1
      nb = nbdcnd
      if (c.eq.0. .and. nb.eq.2) nb = 6
c
c     enter data on theta boundaries
c
      go to (108,108,110,110,112,112,108,110,112),mbdcnd
  108 bm(1) = bm(1)-am(1)
      x = 2.*am(1)
      do 109 j=1,n
         f(1,j) = f(1,j)-x*bda(j)
  109 continue
      go to 112
  110 bm(1) = bm(1)+am(1)
      x = dth*am(1)
      do 111 j=1,n
         f(1,j) = f(1,j)+x*bda(j)
  111 continue
  112 continue
      go to (113,115,115,113,113,115,117,117,117),mbdcnd
  113 bm(m) = bm(m)-cm(m)
      x = 2.*cm(m)
      do 114 j=1,n
         f(m,j) = f(m,j)-x*bdb(j)
  114 continue
      go to 117
  115 bm(m) = bm(m)+cm(m)
      x = dth*cm(m)
      do 116 j=1,n
         f(m,j) = f(m,j)-x*bdb(j)
  116 continue
  117 continue
c
c     enter data on r boundaries
c
      go to (118,118,120,120,122,122),nb
  118 bn(1) = bn(1)-an(1)
      x = 2.*an(1)
      do 119 i=1,m
         f(i,1) = f(i,1)-x*bdc(i)
  119 continue
      go to 122
  120 bn(1) = bn(1)+an(1)
      x = dr*an(1)
      do 121 i=1,m
         f(i,1) = f(i,1)+x*bdc(i)
  121 continue
  122 continue
      go to (123,125,125,123,123,125),nb
  123 bn(n) = bn(n)-cn(n)
      x = 2.*cn(n)
      do 124 i=1,m
         f(i,n) = f(i,n)-x*bdd(i)
  124 continue
      go to 127
  125 bn(n) = bn(n)+cn(n)
      x = dr*cn(n)
      do 126 i=1,m
         f(i,n) = f(i,n)-x*bdd(i)
  126 continue
  127 continue
c
c     check for singular problem.  if singular, perturb f.
c
      pertrb = 0.
      go to (137,137,128,137,137,128,137,128,128),mbdcnd
  128 go to (137,137,129,137,137,129),nb
  129 if (elmbda) 137,131,130
  130 ierr1 = 10
      go to 137
  131 continue
      isw = 2
      do 133 i=1,m
         x = 0.
         do 132 j=1,n
            x = x+f(i,j)
  132    continue
         pertrb = pertrb+x*snth(i)
  133 continue
      x = 0.
      do 134 j=1,n
         x = x+rsq(j)
  134 continue
      pertrb = 2.*(pertrb*sin(dth/2.))/(x*(cos(a)-cos(b)))
      do 136 j=1,n
         x = rsq(j)*pertrb
         do 135 i=1,m
            f(i,j) = f(i,j)-x
  135    continue
  136 continue
  137 continue
      a2 = 0.
      do 138 i=1,m
         a2 = a2+f(i,1)
  138 continue
      a2 = a2/rsq(1)
c
c     initialize blktri
c
      if (intl .ne. 0) go to 139
      call blktri (0,1,n,an,bn,cn,1,m,am,bm,cm,idimf,f,ierr1,wrk)
  139 continue
c
c     call blktri to solve system of equations.
c
      call blktri (1,1,n,an,bn,cn,1,m,am,bm,cm,idimf,f,ierr1,wrk)
      if (isw.ne.2 .or. c.ne.0. .or. nbdcnd.ne.2) go to 143
      a1 = 0.
      a3 = 0.
      do 140 i=1,m
         a1 = a1+snth(i)*f(i,1)
         a3 = a3+snth(i)
  140 continue
      a1 = a1+rsq(1)*a2/2.
      if (mbdcnd .eq. 3)
     1    a1 = a1+(sin(b)*bdb(1)-sin(a)*bda(1))/(2.*(b-a))
      a1 = a1/a3
      a1 = bdc(1)-a1
      do 142 i=1,m
         do 141 j=1,n
            f(i,j) = f(i,j)+a1
  141    continue
  142 continue
  143 continue
      return
      end
