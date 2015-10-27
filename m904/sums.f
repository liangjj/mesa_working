*deck @(#)sums.f	5.1  11/6/94
      subroutine sums(s,t,row)
c.12/01/76
c  access the matrix and compute certain sums for subroutine
c  'simeig'.  some of its details depend on the matrix file
c  structure.  in all cases the matrix is assumed to be arranged by
c  rows of the lower triangle, diagonal elements appearing last in
c  each row, and zero off-diagonal elements omitted (computing time
c  is proportional to the number of elements present).
c
c  there are three normal entries.  entry 'sumsa' is used for the
c  evaluation of column sums during initialization.  entries 'sumsb'
c  and 'sumsc' are used for the evaluation of row sums and column
c  sums, respectively, during the iterations.
c
c  this version brings in one row of the matrix at a time.  each row
c  is preceded by an integer count of the number of words in the row
c  (including the count word), but this count is always read as part
c  of the preceding row.
c
c  each nonzero matrix element a(mu,nu) (mu.ge.nu) occupies one word,
c  with the column index (nu) packed into the least significant bits.
c
c  single-variable arguments are in common /args/, and are described in
c  subroutine 'simeig'.  the use of the array arguments 's' (length=k*n
c  or l*m) and 't' (length=k) depends on the entry.  'row' is used to
c  save the current row of the matrix between 'sumsb' and 'sumsc' calls.
c
      implicit real*8 (a-h,o-z)
      integer and
      common /args/ omega,amumu,bi,n,k,m,l,itmax,it,noconv,mu,idgwrt
      common/tapes/iw,iunt1a,iunt1b,iunt2a,iunts1,iunts2
      common /core/ ecore,thresh,ev,mask,nmask,limit1,limit2,intmxo
*mdc*if cray
*      dimension s(*),t(*),row(*)
*      equivalence (bi,ib)
*mdc*else
      dimension s(*),t(*),row(*),iba(2),nua(2)
      equivalence (iba(1),bi),(iba(2),ib),(nua(1),rnu),(nua(2),nu)
*mdc*endif
c
      entry sumsa(s,t,row)
c
c  this entry adds a(mu,nu)*t(i) to s(i,nu) for nu=1,mu-1 and i=1,k,
c  and places a(mu,mu) in 'amumu'.
c
      if (mu.eq.1) ib=2
      call hin(row,ib)
      nel=ib-2
      bi=row(ib)
      amumu=row(nel+1)
      if (nel.eq.0) return
      if (k-2) 31,32,33
c
 31   do 311 j=1,nel
*mdc*if cray
*        nu= and(row(j),mask)
*mdc*elseif sun
        rnu=row(j)
        nu= and(nu,mask)
*mdc*else
*        rnu=row(j)
*        nu=iand(nu,mask)
*mdc*endif
 311    s(nu)=s(nu)+row(j)*t(1)
      return
c
 32   do 321 j=1,nel
*mdc*if cray
*        nu= and(row(j),mask)
*mdc*elseif sun
        rnu=row(j)
        nu= and(nu,mask)
*mdc*else
*        rnu=row(j)
*        nu=iand(nu,mask)
*mdc*endif
        nu=nu+nu
        s(nu-1)=s(nu-1)+row(j)*t(1)
 321    s(nu)=s(nu)+row(j)*t(2)
      return
c
 33   do 331 j=1,nel
*mdc*if cray
*        nu= and(row(j),mask)
*mdc*elseif sun
        rnu=row(j)
        nu= and(nu,mask)
*mdc*else
*        rnu=row(j)
*        nu=iand(nu,mask)
*mdc*endif
        nu=(nu-1)*k
        do 3311 i=1,k
 3311     s(nu+i)=s(nu+i)+row(j)*t(i)
 331    continue
      return
c
      entry sumsb(s,t,row)
c
c  this entry adds sum(nu=1,mu-1)a(mu,nu)*s(i,nu) to t(i) for i=1,k,
c  and places a(mu,mu) in 'amumu'.
c
      if (mu.eq.1) ib=2
      call hin(row,ib)
      nel=ib-2
      bi=row(ib)
      amumu=row(nel+1)
      if (nel.eq.0) return
      if (k.ne.l) go to 53
      if (k-2) 51,52,53
c
 51   do 511 j=1,nel
*mdc*if cray
*        nu= and(row(j),mask)
*mdc*elseif sun
        rnu=row(j)
        nu= and(nu,mask)
*mdc*else
*        rnu=row(j)
*        nu=iand(nu,mask)
*mdc*endif
 511    t(1)=t(1)+row(j)*s(nu)
      return
c
 52   do 521 j=1,nel
*mdc*if cray
*        nu= and(row(j),mask)
*mdc*elseif sun
        rnu=row(j)
        nu= and(nu,mask)
*mdc*else
*        rnu=row(j)
*        nu=iand(nu,mask)
*mdc*endif
        nu=nu+nu
        t(1)=t(1)+row(j)*s(nu-1)
 521    t(2)=t(2)+row(j)*s(nu)
      return
c
 53   do 531 j=1,nel
*mdc*if cray
*        nu= and(row(j),mask)
*mdc*elseif sun
        rnu=row(j)
        nu= and(nu,mask)
*mdc*else
*        rnu=row(j)
*        nu=iand(nu,mask)
*mdc*endif
        nu=(nu-1)*l
        do 5311 i=1,k
 5311     t(i)=t(i)+row(j)*s(nu+i)
 531    continue
      return
c
      entry sumsc(s,t,row)
c
c  this entry adds a(mu,nu)*t(i) to s(i,nu) for nu=1,mu-1 and i=1,k,
c  using the previously read row of the matrix.
c
      if (k-2) 71,72,73
c
 71   do 711 j=1,nel
*mdc*if cray
*        nu= and(row(j),mask)
*mdc*elseif sun
        rnu=row(j)
        nu= and(nu,mask)
*mdc*else
*        rnu=row(j)
*        nu=iand(nu,mask)
*mdc*endif
 711    s(nu)=s(nu)+row(j)*t(1)
      return
c
 72   do 721 j=1,nel
*mdc*if cray
*        nu= and(row(j),mask)
*mdc*elseif sun
        rnu=row(j)
        nu= and(nu,mask)
*mdc*else
*        rnu=row(j)
*        nu=iand(nu,mask)
*mdc*endif
        nu=nu+nu
        s(nu-1)=s(nu-1)+row(j)*t(1)
 721    s(nu)=s(nu)+row(j)*t(2)
      return
c
 73   do 731 j=1,nel
*mdc*if cray
*        nu= and(row(j),mask)
*mdc*elseif sun
        rnu=row(j)
        nu= and(nu,mask)
*mdc*else
*        rnu=row(j)
*        nu=iand(nu,mask)
*mdc*endif
        nu=(nu-1)*k
        do 7311 i=1,k
 7311     s(nu+i)=s(nu+i)+row(j)*t(i)
 731    continue
      return
c
      end
