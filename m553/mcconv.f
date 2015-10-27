*deck @(#)mcconv.f	5.1  11/6/94
      subroutine mcconv (wo,vec,nts,nrb,wn)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcconv.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8(a-h,o-z)
cc
c     extended dummy wo,vec,wn
cc
      dimension wo(1),vec(1),wn(1)
c
      len=nts*nrb
      do 1 i=1,len
 1    wn(i)=0.d0
c
      k1=0
      l=0
      do 30 i=1,nts
         m=0
         do 20 j=1,nts
            l=l+1
            s1=vec(l)
            k2=k1
            do 20 k=1,nrb
               m=m+1
               k2=k2+1
 20         wn(k2)=wn(k2)+s1*wo(m)
 30      k1=k1+nrb
         return
         end
