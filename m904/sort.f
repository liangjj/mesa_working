*deck @(#)sort.f	5.1  11/6/94
      subroutine sort(a,idx,n)
c
c  sort the array a according to absolute magnitude.
c
c   a         the input array
c  idx        the indexing array
c                  idx(i) is the position of the ith element of the
c                  input array in the ordered array
c
      implicit real*8 (a-h,o-z)
      common/tapes/iw,iunt1a,iunt1b,iunt2a,iunts1,iunts2
      dimension a(*), idx(*)
c
      if(n-1) 90,100,3
c
    3 do 5 i=1,n
    5 idx(i)=i
c
      ns2=n/2
      ns21=ns2+2
      ict=1
      i=2
   10 n1=ns21-i
      nn=n
      ik=n1
   15 c=a(ik)
      ic=idx(ik)
   20 jk=2*ik
      if(jk.gt.nn) go to 35
      if(jk.eq.nn) go to 30
      if(abs(a(jk+1)).ge.abs(a(jk))) go to 30
      jk=jk+1
   30 if(abs(a(jk)).ge.abs(c)) go to 35
      a(ik)=a(jk)
      idx(ik)=idx(jk)
      ik=jk
      go to 20
c
   35 a(ik)=c
      idx(ik)=ic
      if(ict.eq.2) go to 55
      if(i.ge.ns2) go to 45
      i=i+1
      go to 10
c
   45 ict=2
      np2=n+2
      i=2
   50 n1=np2-i
      nn=n1
      ik=1
      go to 15
c
   55 t = a(1)
      a(1) = a(n1)
      a(n1) = t
      it = idx(1)
      idx(1) = idx(n1)
      idx(n1) = it
      if(i.ge.n) return
      i=i+1
      go to 50
c
   90 write (iw,95)
   95 format(///'error return from sort--n less than or equal to 1')
      return
c
  100 idx(1)=1
      return
c
      end
