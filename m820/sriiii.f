*deck @@(#)sriiii.f	5.1  11/6/94
      subroutine sriiii(ni,in,nin,val,lab,bin,lenbin,nbin,
     #     nii,pt,ops,is,bftorb,ijadd,kadd,ladd,
     #     numij,norbs,nsym,isort,rsort,lnsort,unit,calc)
c
      implicit integer (a-z)
c
      character*(*) calc
      character*16 unit
      real*8 in(nii,nin),val(lenbin),rsort(*)
      integer lab(lenbin),bin(lenbin),bftorb(ni),ijadd(numij)
      integer kadd(norbs,nsym),ladd(norbs,nsym),isort(lnsort)
      character*(*) ops
      logical logkey,prnt
c
      common /io/ inp,iout
c
      ioff(i,j)=i*(i-1)/2+j
c
c     -----  check for print option -----
      prnt=logkey(ops,'print=m820=two-electron-integrals',
     $            .false.,' ')
c
c     ----- loop over buffer loads of so integrals -----
c
      k=0
      l=0
      maxkl=0
 1000 continue
         minkl=maxkl+1
         maxkl=min(nii,maxkl+nin)
         numkl=maxkl-minkl+1
         if (calc.eq.'ci') then
            call iosys('read real "mo two-electron integrals" from '//
     $           unit,nii*numkl,in,pt+(minkl-1)*nii,' ')
         else if (calc.eq.'mcscf') then
            call iosys('read real "g mcscf" from '//unit,nii*numkl,in,
     $           pt+(minkl-1)*nii,' ')
         end if
c
c
         if (prnt) then
            write (iout,9001)
 9001       format (5x,'m820: two-electron integrals (mo)')
            call matout(in,nii,numkl,nii,numkl,iout)
         end if
c
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
            korb=bftorb(k)
            lorb=bftorb(l)
c
            ij=0
            do 2 i=1,ni
               iorb=bftorb(i)
               if (nbin+i.gt.lenbin) then
                  call sorter('with bin',isort,rsort,0,nbin,lab,
     #                         bin,val,0,0,0,0)
                  nbin=0
               end if
               do 1 j=1,i
                  ij=ij+1
                  if (ij.lt.minkl+kl-1) go to 1
                  jorb=bftorb(j)
                  nbin=nbin+1
                  lab(nbin)=guglab(iorb,is,jorb,is,korb,is,lorb,is,
     #                             ijadd,numij,kadd,ladd,norbs,nsym)
                  val(nbin)=in(ij,kl)
    1          continue
    2       continue
   10    continue
      if (maxkl.lt.nii) go to 1000
c
c
      return
      end
