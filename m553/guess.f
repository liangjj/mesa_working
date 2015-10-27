*deck @(#)guess.f	5.1  11/6/94
      subroutine guess(buf,diag,smlham,ipt,n,nsml,nlook)
c
c***begin prologue     guess
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)guess.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       guess
c
      implicit real*8(a-h,o-z)
c
      dimension buf(n,n),smlham(2),diag(2),ipt(2)
c
      common /io/ inp,iout
c
      big=1.d+25
      small=1.d+20
      do 20 i=1,nsml
         xmin=small
         imin=0
         do 10 j=1,nlook
            if(diag(j).gt.xmin) go to 10
            imin=j
            xmin=diag(j)
 10      continue
         if(imin.eq.0)go to 1000
         ipt(i)=imin
         diag(imin)=big
 20   continue
c
      ix=0
      do 120 i=1,nsml
         ii=ipt(i)
         diag(ii)=buf(ii,ii)
         do 110 j=1,i
            ix=ix+1
            jj=ipt(j)
            smlham(ix)=buf(ii,jj)
 110     continue
 120  continue
cc
      return
cc
 1000 continue
      write(iout,1010)
 1010 format('  imin = 0 in guess   bug .. stop ')
      call lnkerr(' ')
c
c
      end
