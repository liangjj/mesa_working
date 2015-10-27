*deck @(#)intgrl.f	5.1  11/6/94
      subroutine intgrl(bfsym,orbtbf,kadd,ladd,ijadd,ijgrp,inint
     #,                 inext,jmnnxt,jmxnxt,ningrp,ijxx,klxx,nklxx
     #,                 ijww,klww,nklww,no34x)
c
c***********************************************************************
c                                                                      *
c pws 19 july 1982                                                     *
c     modified 19 july 1982 to form 3 and 4-external addressing        *
c     scheme and to generate ij-blocks rather than only i-blocks       *
c                                                                      *
c     count and compute addresses for the integrals including symmetry.*
c     the address of an integral with indices i>j>k>l (or equal) is    *
c     given by the following expression:                               *
c                                                                      *
c ijadd(ij)+kadd(k+sym(i)sym(l)*norbs)+ladd(l+sym(i)sym(l)sym(k)*norbs)*
c                                                                      *
c     where ij=i*(i-1)/2+j and norbs is number of orbitals in ci. note *
c     also that symmetries need direct products. the integrals are     *
c     stored with the following offset from the address above:         *
c                                                                      *
c       type            1            2          3                      *
c         1          (ik,jl)      (ij,kl)    (il,jk)                   *
c         2          (ij,jl)      (il,jj)                              *
c         3          (ik,il)      (ii,kl)                              *
c         4          (il,jl)      (ij,ll)                              *
c         5          (ii,il)      (il,ll)    <i/h/l>                   *
c         6          (il,il)      (ii,ll)                              *
c         7          (ii,ii)      <i/h/i>                              *
c                                                                      *
c     the address of 3- and 4-external matrix elements is given by     *
c                                                                      *
c             wy entry                   xy entry                      *
c        ijww(ia)+klww(bc)    or    ijxx(ia)+klxx(bc)       (3x)       *
c                             or                                       *
c        ijww(ab)+klww(cd)    or    ijxx(ab)+klxx(cd)       (4x)       *
c             ww entry                   xx entry                      *
c                                                                      *
c     this subroutine also determines the integral block size (nmax).  *
c     ijgrp(ij) gives block containing all integrals with i and j      *
c     indices.                                                         *
c***********************************************************************
c
c
      implicit integer (a-z)
      integer numint
      logical no34x
c
      common /bloksz/ blksiz,absmax,maxsiz
      common /tapes/  out,errout,input,drttap
      common /dimens/ nbf,nsym,norbs,nrowsp,nrows4p,nrows,nrows4
     #,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     #,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc,nvref
     #,               nijvir
c
      dimension bfsym(nbf),orbtbf(norbs),kadd(symorb),ladd(symorb)
      dimension ijadd(numij),ijgrp(numij),inint(norbs)
      dimension inext(norbs),jmnnxt(norbs),jmxnxt(norbs),ningrp(1)
      dimension ijxx(numij),klxx(numij),ijww(numij)
      dimension klww(numij),nklww(nsym,orbfrm),nklxx(nsym,orbfrm)
c
   21 format(5x,'integral storage:',
     $      /8x,'number of groups         ',i10,
     $      /8x,'group size               ',i10,
     $      /8x,'number of integrals      ',i10)
c
c
      do 101 i=1,numij
         ijxx(i)=0
         ijww(i)=0
  101 continue
c
c     ----- count integrals and form kadd and ladd arrays -----
c
      numint=0
      jadmax=0
      ijad=0
      do 9 i=1,norbs
         isym=bfsym(orbtbf(i))
         jad=0
         ij=i*(i-1)/2
         do 8 j=1,i
            ij=ij+1
            ijsym=xor(isym,bfsym(orbtbf(j)))
            kad=0
            do 7 k=1,j
               ijkad=k+ijsym*norbs
               if (kadd(ijkad).eq.0.or.kadd(ijkad).eq.kad) go to 2
               write (errout,1)
    1          format (//,' symmetry problems with k in intgrl',//)
               call lnkerr('symmetry problems with k index of '//
     #              'integrals')
    2          continue
               kadd(ijkad)=kad
               ijksym=xor(ijsym,bfsym(orbtbf(k)))
               lad=0
               do 6 l=1,k
                  ijklad=l+ijksym*norbs
                  if (ladd(ijklad).eq.0.or.ladd(ijklad).eq.lad) go to 4
                  write (errout,3)
    3             format (//,' symmetry problems with l in intgrl',//)
                  call lnkerr('symmetry problems with l index of'//
     #                        ' integrals')
    4             continue
                  ladd(ijklad)=lad
                  if (ijksym.ne.bfsym(orbtbf(l))) go to 5
                  if (l.eq.j.and.l.lt.i) go to 5
                  lad=lad+3
c     if (k.eq.l) lad=lad-1
    5             continue
    6          continue
               kad=kad+lad
    7       continue
            jad=jad+kad
            ijadd(ij)=kad
    8    continue
         numint=numint+jad
         inint(i)=jad
         if (jad.gt.jadmax) jadmax=jad
    9 continue
      if(no34x) goto 300
c
c     ----- work out kl addressing for 3- and 4- external integrals ----
c
      do 34 isym=1,nsym
         num=0
         do 33 i=2,orbfrm
            is=xor(bfsym(orbtbf(i)),(isym-1))
            ia=i*(i-1)/2
            do 32 j=1,i-1
               if (bfsym(orbtbf(j)).ne.is) go to 31
               num=num+1
               klxx(ia+j)=num
   31          continue
   32       continue
            nklxx(isym,i)=num
   33    continue
c     nklxx(isym)=num
   34 continue
c
      do 38 isym=1,nsym
         num=0
         do 37 i=1,orbfrm
            is=xor((isym-1),bfsym(orbtbf(i)))
            ia=i*(i-1)/2
            do 36 j=1,i
               if (bfsym(orbtbf(j)).ne.is) go to 35
               num=num+1
               klww(ia+j)=num
   35          continue
   36       continue
            nklww(isym,i)=num
   37    continue
c     nklww(isym)=num
   38 continue
c
c     ----- number of 3- and 4-external elements in ij blocks -----
c
      do 44 junk=orbfrm+1,norbs
         i=norbs-junk+orbfrm+1
         ia=i*(i-1)/2
         is=bfsym(orbtbf(i))
         num=0
         do 42 j=1,orbfrm
            ijs=xor(is,bfsym(orbtbf(j)))
            ijww(ia+j)=nklww(ijs+1,orbfrm)
            num=num+nklww(ijs+1,orbfrm)
   42    continue
         do 43 j=1,orbfrm
            ijs=xor(is,bfsym(orbtbf(j)))
            ijxx(ia+j)=nklxx(ijs+1,orbfrm)
            num=num+nklxx(ijs+1,orbfrm)
   43    continue
         inint(i)=inint(i)+num
         numint=numint+num
   44 continue
c
      do 41 junk=1,orbfrm
         i=orbfrm-junk+1
         num=0
         ia=i*(i-1)/2
         is=bfsym(orbtbf(i))
         do 39 j=1,i-1
            ijs=xor(is,bfsym(orbtbf(j)))
            ijxx(ia+j)=nklxx(ijs+1,orbfrm)
            num=num+nklxx(ijs+1,orbfrm)
   39    continue
c
         do 40 j=1,i
            ijs=xor(is,bfsym(orbtbf(j)))
            ijww(ia+j)=nklww(ijs+1,orbfrm)
            num=num+nklww(ijs+1,orbfrm)
   40    continue
         inint(i)=inint(i)+num
         numint=numint+num
   41 continue
c
c     ----- work out a convenient output block size for integrals -----
c
  300 continue
      jadmax=0
      do 45 i=1,norbs
         if (inint(i).gt.jadmax) jadmax=inint(i)
   45 continue
      nmax=blksiz
      if (maxsiz-2*nwks.lt.blksiz) nmax=maxsiz-2*nwks
      if (nmax.lt.jadmax) nmax=blksiz
      if(nmax.le.jadmax) then
         write(out,10) jadmax
   10    format(5x,'necessary i-block size:',8x,i7,
     $         /8x,'will try to fit in largest i-block.')
c  10 format (//,' cannot fit i-block of integrals in requested block'
c    #,          ' size',/,' need',i7,' integrals per block, so will'
c    #,          ' try to hold largest i-block',//)
         nmax=jadmax
      endif
      if(nmax.gt.absmax) then
         write(out,12) nmax, absmax
   12    format(8x,'cannot handle largest i-block. size: ',i8,
     $         /8x,'maximum block size                 : ',i8,
     $         /8x,'will split i-blocks.')
         nmax=absmax
      endif
c  12 format (//,' cannot handle i-block of integrals of size',i8
c    #,        /,' so will split i-blocks',//)
      if (nmax.gt.numint) nmax=numint
c
c     ----- generate the ijadd and ijgrp arrays -----
c
      group=0
      left=0
      sum=0
      do 59 i=norbs,1,-1
         ia=i*(i-1)/2
         if (inint(i).le.left) go to 50
         if (group.gt.0) ningrp(group)=sum
         group=group+1
         left=nmax
         sum=0
   50    continue
         ijmax=ia+i
         do 55 j=i,1,-1
            ij=ia+j
            nij=ijadd(ij)+ijww(ij)+ijxx(ij)
            if(no34x) nij=ijadd(ij)
            if (nij.gt.nmax) go to 910
            if (nij.le.left) go to 54
            ijmin=ij+1
            do 51 ijq=ijmax,ijmin,-1
               t=ijadd(ijq)
               ijadd(ijq)=sum
               sum=sum+t
   51       continue
            if(.not.no34x) then
               do 52 ijq=ijmax,ijmin,-1
                  t=ijww(ijq)
                  ijww(ijq)=sum
                  sum=sum+t
   52          continue
               do 53 ijq=ijmax,ijmin,-1
                  t=ijxx(ijq)
                  ijxx(ijq)=sum
                  sum=sum+t
   53          continue
            end if
c
            ijmax=ij
            ningrp(group)=sum
            group=group+1
            sum=0
            left=nmax
   54       continue
            left=left-nij
            ijgrp(ij)=group
   55    continue
c
         ijmin=ij
         do 56 ijq=ijmax,ijmin,-1
            t=ijadd(ijq)
            ijadd(ijq)=sum
            sum=sum+t
   56    continue
         if(.not.no34x) then
            do 57 ijq=ijmax,ijmin,-1
               t=ijww(ijq)
               ijww(ijq)=sum
               sum=sum+t
   57       continue
            do 58 ijq=ijmax,ijmin,-1
               t=ijxx(ijq)
               ijxx(ijq)=sum
               sum=sum+t
   58       continue
         end if
c
   59 continue
      ngroup=group
      write (out,21) ngroup,nmax,numint
c
c     ----- generate inext and jnext arrays -----
c
      do 20 orb=1,norbs
         next=norbs-orb+1
         inext(next)=orb
         jmnnxt(next)=1
         jmxnxt(next)=orb
   20 continue
      return
c
  910 continue
      write (errout,911) nmax,nij
  911 format (//,' cannot fit ij-block of integrals. nmax=',i7
     #,        / '                        size of ij-block',i7,//)
      call lnkerr('cannot fit ij-block of integrals')
      end
