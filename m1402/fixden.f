*deck @(#)fixden.f	5.1  11/6/94
      subroutine fixden(den,guga1,orbtbf,nbf)
      implicit integer(a-z)
      common /io/ inp,iout
c
      real*8 den(nbf,nbf),guga1(nbf,nbf)
      integer orbtbf(nbf)
c
      do 1 i=1,nbf
         ic=orbtbf(i)
         do 2 j=1,nbf
            jc=orbtbf(j)
            den(ic,jc)=guga1(i,j)
 2       continue
 1    continue
c
      return
      end
