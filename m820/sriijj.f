*deck @@(#)sriijj.f	5.1  11/6/94
      subroutine sriijj(ni,nj,in,nin,val,lab,bin,lenbin,nbin,
     #                  nii,njj,pt,ops,is,js,ibforb,jbforb,
     #                  ijadd,kadd,ladd,numij,norbs,nsym,
     #                  isort,rsort,lnsort,unit,calc)
c
      implicit integer (a-z)
c
      character*(*) calc
      character*16 unit
      real*8 in(nii,nin),val(lenbin),rsort(*)
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
      l=0
      maxkl=0
 1000 continue
         minkl=maxkl+1
         maxkl=min(njj,maxkl+nin)
         numkl=maxkl-minkl+1
         if (calc.eq.'ci') then
            call iosys('read real "mo two-electron integrals" from '//
     $           unit,nii*numkl,in,pt+(minkl-1)*nii,' ')
         else if (calc.eq.'mcscf') then
            call iosys('read real "g mcscf" from '//unit,nii*numkl,in,
     $           pt+(minkl-1)*nii,' ')
         end if
c
c        ----- sift through a canonical list -----
c
         do 10 kl=1,numkl
c
            l=l+1
            if (l.gt.k) then
               k=k+1
               l=1
            end if
            korb=jbforb(k)
            lorb=jbforb(l)
c
            ij=0
            do 2 i=1,ni
               iorb=ibforb(i)
               if (nbin+i.gt.lenbin) then
                  call sorter('with bin',isort,rsort,0,nbin,lab,
     #                         bin,val,0,0,0,0)
                  nbin=0
               end if
               do 1 j=1,i
                  ij=ij+1
                  jorb=ibforb(j)
                  nbin=nbin+1
                  lab(nbin)=guglab(iorb,is,jorb,is,korb,js,lorb,js,
     #                             ijadd,numij,kadd,ladd,norbs,nsym)
                  val(nbin)=in(ij,kl)
    1          continue
    2       continue
   10    continue
      if (maxkl.lt.njj) go to 1000
c
c
      return
      end
