*deck @(#)mcstvt.f	5.1  11/6/94
      subroutine mcstvt(result,vector,scalar,nr,ns)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcstvt.f	5.1   11/6/94
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
cmp   extended dummy result,vector
cc
c
c--------------------------------------------------------------------c
c    the transpose of the vector array is used in this dot product
c--------------------------------------------------------------------c
c
      dimension result(2),vector(2)
c
      nxe=0
      do 20 j=1,ns
         jx=j
         nxs=nxe+1
         nxe=nxe+nr
         do 10 i=nxs,nxe
            result(i)=result(i)+scalar*vector(jx)
            jx=jx+nr
 10      continue
 20   continue
c
      return
      end
