*deck poisp2
      subroutine poisp2 (m, n, a, bb, c, q, idimq, b, b2, b3, w, w2, w3,
     +   d, tcos, p)
c***begin prologue  poisp2
c***subsidiary
c***purpose  subsidiary to genbun
c***library   slatec
c***type      single precision (poisp2-s, cmposp-c)
c***author  (unknown)
c***description
c
c     subroutine to solve poisson equation with periodic boundary
c     conditions.
c
c***see also  genbun
c***routines called  poisd2, poisn2
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  poisp2
c
      dimension       a(*)       ,bb(*)      ,c(*)       ,q(idimq,*) ,
     1                b(*)       ,b2(*)      ,b3(*)      ,w(*)       ,
     2                w2(*)      ,w3(*)      ,d(*)       ,tcos(*)    ,
     3                p(*)
c***first executable statement  poisp2
      mr = m
      nr = (n+1)/2
      nrm1 = nr-1
      if (2*nr .ne. n) go to 107
c
c     even number of unknowns
c
      do 102 j=1,nrm1
         nrmj = nr-j
         nrpj = nr+j
         do 101 i=1,mr
            s = q(i,nrmj)-q(i,nrpj)
            t = q(i,nrmj)+q(i,nrpj)
            q(i,nrmj) = s
            q(i,nrpj) = t
  101    continue
  102 continue
      do 103 i=1,mr
         q(i,nr) = 2.*q(i,nr)
         q(i,n) = 2.*q(i,n)
  103 continue
      call poisd2 (mr,nrm1,1,a,bb,c,q,idimq,b,w,d,tcos,p)
      ipstor = w(1)
      call poisn2 (mr,nr+1,1,1,a,bb,c,q(1,nr),idimq,b,b2,b3,w,w2,w3,d,
     1             tcos,p)
      ipstor = max(ipstor,int(w(1)))
      do 105 j=1,nrm1
         nrmj = nr-j
         nrpj = nr+j
         do 104 i=1,mr
            s = .5*(q(i,nrpj)+q(i,nrmj))
            t = .5*(q(i,nrpj)-q(i,nrmj))
            q(i,nrmj) = s
            q(i,nrpj) = t
  104    continue
  105 continue
      do 106 i=1,mr
         q(i,nr) = .5*q(i,nr)
         q(i,n) = .5*q(i,n)
  106 continue
      go to 118
  107 continue
c
c     odd  number of unknowns
c
      do 109 j=1,nrm1
         nrpj = n+1-j
         do 108 i=1,mr
            s = q(i,j)-q(i,nrpj)
            t = q(i,j)+q(i,nrpj)
            q(i,j) = s
            q(i,nrpj) = t
  108    continue
  109 continue
      do 110 i=1,mr
         q(i,nr) = 2.*q(i,nr)
  110 continue
      lh = nrm1/2
      do 112 j=1,lh
         nrmj = nr-j
         do 111 i=1,mr
            s = q(i,j)
            q(i,j) = q(i,nrmj)
            q(i,nrmj) = s
  111    continue
  112 continue
      call poisd2 (mr,nrm1,2,a,bb,c,q,idimq,b,w,d,tcos,p)
      ipstor = w(1)
      call poisn2 (mr,nr,2,1,a,bb,c,q(1,nr),idimq,b,b2,b3,w,w2,w3,d,
     1             tcos,p)
      ipstor = max(ipstor,int(w(1)))
      do 114 j=1,nrm1
         nrpj = nr+j
         do 113 i=1,mr
            s = .5*(q(i,nrpj)+q(i,j))
            t = .5*(q(i,nrpj)-q(i,j))
            q(i,nrpj) = t
            q(i,j) = s
  113    continue
  114 continue
      do 115 i=1,mr
         q(i,nr) = .5*q(i,nr)
  115 continue
      do 117 j=1,lh
         nrmj = nr-j
         do 116 i=1,mr
            s = q(i,j)
            q(i,j) = q(i,nrmj)
            q(i,nrmj) = s
  116    continue
  117 continue
  118 continue
c
c     return storage requirements for p vectors.
c
      w(1) = ipstor
      return
      end
