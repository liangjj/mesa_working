*deck @(#)fixval.f	1.1  11/30/90
      subroutine fixval(values,xval,nao)
c
c***begin prologue     fixval
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)fixval.f	1.1   11/30/90
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       fixval
c
      implicit real*8 (a-h,o-z)
      dimension values(2),xval(2)
c
      nao2=(nao*(nao+1))/2
      nao4=nao*nao
      ix=0
      ii=0
      do 20 i=1,nao
         do 10 j=1,i
            ix=ix+1
            ij=ii+j
            call scopy(nao2,xval(ij),nao4,values(ix),nao2)
 10      continue
         ii=ii+nao
 20   continue
c
      return
      end
