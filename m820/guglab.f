*deck @@(#)guglab.f	5.1  11/6/94
      function guglab(it,ist,jt,jst,kt,kst,lt,lst,ijadd,numij,
     #                        kadd,ladd,norbs,nsym)
c
      implicit integer (a-z)
      integer guglab
c
      integer ijadd(numij),kadd(norbs,nsym),ladd(norbs,nsym)
      integer temp(4,2)
c
      equivalence (temp(1,1),i),(temp(2,1),j)
      equivalence (temp(3,1),k),(temp(4,1),l)
      equivalence (temp(1,2),is),(temp(2,2),js)
      equivalence (temp(3,2),ks),(temp(4,2),ls)
c
      ioff(i,j)=i*(i-1)/2+j
c
c     ----- make sure we have a canonical representation -----
c
      ic=it
      jc=jt
      kc=kt
      lc=lt
c
      if (jc.gt.ic) then
         t=ic
         ic=jc
         jc=t
      end if
      if (lc.gt.kc) then
         t=kc
         kc=lc
         lc=t
      end if
      if (kc.gt.ic.or.(ic.eq.kc.and.lc.gt.jc)) then
         t=ic
         ic=kc
         kc=t
         t=jc
         jc=lc
         lc=t
      end if
c
      temp(1,1)=it
      temp(2,1)=jt
      temp(3,1)=kt
      temp(4,1)=lt
      temp(1,2)=ist
      temp(2,2)=jst
      temp(3,2)=kst
      temp(4,2)=lst
c
      pass=0
    1 continue
         diff=0
         do 3 n=2,4
            if (temp(n,1).gt.temp(n-1,1)) then
               t=temp(n,1)
               temp(n,1)=temp(n-1,1)
               temp(n-1,1)=t
               t=temp(n,2)
               temp(n,2)=temp(n-1,2)
               temp(n-1,2)=t
               diff=diff+1
            end if
    3    continue
         pass=pass+1
         if (pass.gt.50) call lnkerr('bubble sort')
      if (diff.gt.0) go to 1
c
      ijs=xor(is-1,js-1)+1
      ijks=xor(ijs-1,ks-1)+1
      if (ijks.ne.ls) then
         guglab=-9999999
         go to 1000
      end if
      ij=ioff(i,j)
c
      guglab=ijadd(ij)+kadd(k,ijs)+ladd(l,ijks)+1
      if (i.eq.l) then
         go to 1000
      else if (i.eq.k) then
         go to 1000
      else if (i.eq.j) then
         if (ic.eq.jc) guglab=guglab+1
         go to 1000
      else
         if (j.eq.l) then
c
c          ----- special case: [il;ll] -----
c
            ij=ioff(i,i)
            if (is.ne.js) then
               guglab=-9999998
               go to 1000
            end if
            guglab=ijadd(ij)+kadd(i,1)+ladd(l,is)+2
         else if (j.eq.k) then
            if (kc.eq.lc) guglab=guglab+1
            go to 1000
         else if (k.eq.l) then
            if (kc.eq.lc) guglab=guglab+1
            go to 1000
         else
            if (j.eq.jc) then
               guglab=guglab+1
               go to 1000
            else if (l.eq.jc) then
               guglab=guglab+2
               go to 1000
            end if
         end if
         go to 1000
      end if
c
c
 1000 continue
c      write (0,1001) ic,jc,kc,lc,it,jt,kt,lt,ist,jst,kst,lst,
c     #               i,j,k,l,is,js,ks,ls,guglab
c 1001 format (1x,4i2,2x,4i2,2x,4i1,5x,4i2,2x,4i1,5x,i10)
c
c
      return
      end
