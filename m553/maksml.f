*deck @(#)maksml.f	5.1  11/6/94
      subroutine maksml(smlham,c,b,n,mdim,nvec)
c
c***begin prologue     maksml
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)maksml.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       maksml
c
      implicit real*8(a-h,o-z)
c
      dimension smlham(2),c(2),b(2)
c
      mold=mdim-nvec
      ix=(mold*(mold+1))/2
      ioff=mold*n
      do 30 m=1,nvec
         js=ioff+1
         je=ioff+n
         ij=0
         iend=mold+m
         do 20 i=1,iend
            xx=0.d0
            do 10 j=js,je
               ij=ij+1
               xx=xx+c(ij)*b(j)
 10         continue
            ix=ix+1
            smlham(ix)=xx
 20      continue
         ioff=ioff+n
 30   continue
c
      return
      end
