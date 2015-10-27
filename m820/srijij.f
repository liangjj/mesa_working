*deck @@(#)srijij.f	5.1  11/6/94
      subroutine srijij(ni,nj,in,nin,val,lab,bin,lenbin,nbin,
     #                  nij,pt,ops,is,js,ibforb,jbforb,
     #                  ijadd,kadd,ladd,numij,norbs,nsym,
     #                  isort,rsort,lnsort,unit,calc)
c
      implicit integer (a-z)
c
      character*(*) calc
      character*16 unit
      real*8 in(nij,nin),val(lenbin),rsort(*)
      integer isort(lnsort)
      integer lab(lenbin),bin(lenbin),ibforb(ni),jbforb(nj)
      integer ijadd(numij),kadd(norbs,nsym),ladd(norbs,nsym)
      character*(*) ops
c
      ioff(i,j)=i*(i-1)/2+j
c
c     ----- loop over buffer loads of so integrals -----
c
      k=0
      l=1
      maxkl=0
 1000 continue
         minkl=maxkl+1
         maxkl=min(nij,maxkl+nin)
         numkl=maxkl-minkl+1
         if (calc.eq.'ci') then
            call iosys('read real "mo two-electron integrals" from '//
     $           unit,nij*numkl,in,pt+(minkl-1)*nij,' ')
         else if (calc.eq.'mcscf') then
            call iosys('read real "g mcscf" from '//unit,nij*numkl,in,
     $           pt+(minkl-1)*nij,' ')
         end if
c
c        ----- sift through a canonical list -----
c
         do 10 kl=1,numkl
c
            k=k+1
            if (k.gt.ni) then
               l=l+1
               k=1
            end if
            korb=ibforb(k)
            lorb=jbforb(l)
c
            ij=0
            do 2 j=1,nj
               jorb=jbforb(j)
               if (nbin+ni.gt.lenbin) then
                  call sorter('with bin',isort,rsort,0,nbin,lab,
     #                         bin,val,0,0,0,0)
                  nbin=0
               end if
               do 1 i=1,ni
                  ij=ij+1
                  if (ij.lt.minkl+kl-1) go to 1
                  iorb=ibforb(i)
                  nbin=nbin+1
                  lab(nbin)=guglab(iorb,is,jorb,js,korb,is,lorb,js,
     #                             ijadd,numij,kadd,ladd,norbs,nsym)
                  val(nbin)=in(ij,kl)
    1          continue
    2       continue
   10    continue
      if (maxkl.lt.nij) go to 1000
c
c
      return
      end
