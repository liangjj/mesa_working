*deck @(#)qntgrl.f	1.2  7/30/91
      subroutine qntgrl(norbs,nsym,nnp,orbsym,ijgrp,ijadd,kadd,ladd,
     #                  inint,ijxx,ijww,klxx,klww,nklxx,nklww,
     #                  nnpvir,orbfrm,out,nmax,ngroup)
c
c***begin prologue  intgrl
c***date written   850109   (yymmdd)
c***revision date  yymmdd   (yymmdd)
c***keywords  distinct row table, drt, integral storage
c
c***author  saxe, paul,    (lanl)
c***purpose  to determine the integral-addressing arrays for the ci
c
c***description
c
c
c***references
c
c***routines called  (none)
c***end prologue  intgrl
c
      implicit integer (a-z)
c
      integer orbsym(norbs),ijgrp(nnp),ijadd(nnp),kadd(norbs,nsym)
      integer ladd(norbs,nsym),inint(norbs)
      integer ijxx(nnp),ijww(nnp),klxx(nnpvir),klww(nnpvir)
      integer nklxx(nsym,orbfrm),nklww(nsym,orbfrm)
c
c     ----- form the offsets for the last (or l) index -----
c
      do 2 sym=1,nsym
         lad=0
         do 1 orb=1,norbs
            if (sym.eq.orbsym(orb)) then
               ladd(orb,sym)=lad
               lad=lad+3
            else
               ladd(orb,sym)=-7999999
            end if
    1    continue
    2 continue
c
c     ----- form offsets for the third (or k) index -----
c
      do 4 sym=1,nsym
         kad=0
         do 3 orb=1,norbs
            kad=kad+max(0,ladd(orb,xor(sym-1,orbsym(orb)-1)+1))
            kadd(orb,sym)=kad
    3    continue
    4 continue
c
c     ----- and the location of the ij blocks -----
c
      ij=0
      do 9 i=1,norbs
         isym=orbsym(i)
         inint(i)=0
         do 8 j=1,i
            ijad=0
            ij=ij+1
            ijsym=xor(isym-1,orbsym(j)-1)+1
            do 6 k=j,1,-1
               ijksym=xor(ijsym-1,orbsym(k)-1)+1
               do 5 l=k,1,-1
                  if (orbsym(l).eq.ijksym.and.
     #                 (i.eq.j.or.l.ne.j)) go to 7
    5          continue
    6       continue
            go to 20
    7       continue
            ijad=ijad+kadd(k,ijsym)+ladd(l,ijksym)+3
   20       continue
            ijadd(ij)=ijad
            inint(i)=inint(i)+ijad
    8    continue
    9 continue
c
c     ----- work out kl addressing for 3- and 4-external supermatrices
c
      do 25 isym=1,nsym
         num=0
         do 24 i=2,orbfrm
            is=xor(orbsym(i)-1,isym-1)+1
            ia=i*(i-1)/2
            do 23 j=1,i-1
               if (orbsym(j).eq.is) then
                  num=num+1
                  klxx(ia+j)=num
               end if
   23       continue
            nklxx(isym,i)=num
   24    continue
   25 continue
c
      do 35 isym=1,nsym
         num=0
         do 34 i=1,orbfrm
            is=xor(orbsym(i)-1,isym-1)+1
            ia=i*(i-1)/2
            do 33 j=1,i
               if (orbsym(j).eq.is) then
                  num=num+1
                  klww(ia+j)=num
               end if
   33       continue
            nklww(isym,i)=num
   34    continue
   35 continue
c
c     ----- number of 3- and 4-external elements per ij block
c
      call izero(ijxx,nnp)
      call izero(ijww,nnp)
c
      if (orbfrm.le.0) go to 503
c
      do 43 i=norbs,1,-1
         ia=i*(i-1)/2
         is=orbsym(i)-1
         num=0
         do 40 j=1,min(orbfrm,i)
            ijs=xor(is,orbsym(j)-1)+1
            ijww(ia+j)=nklww(ijs,orbfrm)
            num=num+nklww(ijs,orbfrm)
   40    continue
         do 41 j=1,min(orbfrm,i-1)
            ijs=xor(is,orbsym(j)-1)+1
            ijxx(ia+j)=nklxx(ijs,orbfrm)
            num=num+nklxx(ijs,orbfrm)
   41    continue
         inint(i)=inint(i)+num
   43 continue
  503 continue
c
c     ----- now we know how many integrals are in each section,
c           form total offsets
c
      sum=0
      ij=(norbs+1)*norbs/2+1
      do 50 i=norbs,1,-1
         ijsv=ij
         do 47 j=i,1,-1
            ij=ij-1
            t=max(0,ijadd(ij))
            ijadd(ij)=sum
            sum=sum+t
   47    continue
            if (orbfrm.le.0) go to 504
         ij=ijsv
         do 48 j=i,1,-1
            ij=ij-1
            t=ijww(ij)
            ijww(ij)=sum
            sum=sum+t
   48    continue
         ij=ijsv
         do 49 j=i,1,-1
            ij=ij-1
            t=ijxx(ij)
            ijxx(ij)=sum
            sum=sum+t
   49    continue
  504    continue
   50 continue
ctemp
      do 51 i=1,nnp
         ijgrp(i)=1
   51 continue
      nmax=sum
      ngroup=1
cend
c
c
      return
      end
